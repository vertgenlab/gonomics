package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

type Settings struct {
	InFile           string
	OutFile          string
	AlignQualFilter  int
	AlignLenFilter   int
	FilterCigar      string
	CollapseByUmi    bool
	SingleCellFormat bool
	FilterByRegion   string
	FilterByFlag     int
}

func filterByQuality(a sam.Sam, filter int) bool {
	if int(a.MapQ) > filter {
		return true
	} else {
		return false
	}
}

func filterByLength(a sam.Sam, filter int) bool {
	var pass bool = false
	for _, k := range a.Cigar {
		if k.Op == 77 && k.RunLength >= filter {
			pass = true
		}
	}
	return pass
}

func filterByCigar(a sam.Sam, filter string) bool {
	if filter == "starrSeqIntrons" {
		for _, k := range a.Cigar { //filter out introns longer than 500bp (the length of 1 standard STARR-seq construct)
			if k.Op == 78 && k.RunLength > 500 {
				return false
			}
		}
	} else {
		if cigar.ToString(a.Cigar) != filter {
			return false
		}
	}
	return true
}

func filterByFlag(a sam.Sam, filter int) bool {
	if int(a.Flag) != filter {
		return false
	}
	return true
}

func filterByRegion(a sam.Sam, filter bed.Bed) bool {
	if bed.Overlap(filter, convert.SamToBed(a)) {
		return true
	} else {
		return false
	}
}

func collapseUMI(a sam.Sam) map[string]sam.Sam {
	var umiMap = make(map[string]sam.Sam)
	sc := sam.ToSingleCellAlignment(a)
	umiBxString := fmt.Sprintf("%s%s", dna.BasesToString(sc.Bx), dna.BasesToString(sc.Umi))
	_, found := umiMap[umiBxString]
	if !found { //if the UMI isn't in the map, add it.
		umiMap[umiBxString] = a
	} else {
		if umiMap[umiBxString].MapQ > a.MapQ {
			//if the UMI is in the map and the UMI in the map has a better alignment score -- do nothing.
		} else if umiMap[umiBxString].MapQ == a.MapQ {
			if umiMap[umiBxString].Pos >= a.Pos { //if the UMIs are tied on alignment score, prioritize the one with a higher start position since the construct barcodes are on the end of constructs,
				//higher chance that the read contains the construct barcode.
				//do nothing
			} else if umiMap[umiBxString].Pos < a.Pos {
				umiMap[umiBxString] = a
			}
		} else if umiMap[umiBxString].MapQ < a.MapQ { //if the UMI in the map has a lower alignment score than the query, replace it.
			umiMap[umiBxString] = a
		}
	}
	return umiMap
}

func appendUmiCbc(a sam.Sam) (sam.Sam, bool) {
	var found bool = false
	var newReadName string

	cbc, cbcFound, _ := sam.QueryTag(a, "CB")
	umi, umiFound, _ := sam.QueryTag(a, "UB")
	cbcString := fmt.Sprintf("%s", cbc)
	w := strings.Split(cbcString, "-")

	if cbcFound && umiFound {
		found = true
		newReadName = fmt.Sprintf("%s_UMI:%s_BX:%s", a.QName, umi, w[0]) //compatible with scCount.go
		a.QName = newReadName
	}
	return a, found
}

func runFilter(s Settings) {
	var passQual, passLen, passCigar, foundUMI, passFlag, passLocation bool
	var newSamMap = make(map[string]sam.Sam)

	samChan, header := sam.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)
	sam.WriteHeaderToFileHandle(out, header)

	var words, p []string
	var chrom string
	var start, end int
	var bedRegion bed.Bed

	if s.FilterByRegion != "" {
		words = strings.Split(s.FilterByRegion, ":")
		chrom = words[0]
		p = strings.Split(words[1], "-")
		if len(p) == 2 {
			start = common.StringToInt(p[0])
			end = common.StringToInt(p[1])
			bedRegion = bed.Bed{
				Chrom:      chrom,
				ChromStart: start,
				ChromEnd:   end,
			}
		} else {
			bedRegion = bed.Bed{
				Chrom:      chrom,
				ChromStart: 0,
				ChromEnd:   numbers.MaxInt,
			}
		}

		for i := range samChan {
			if s.AlignQualFilter > 0 {
				passQual = filterByQuality(i, s.AlignQualFilter)
				if !passQual {
					continue
				}
			}
			if s.AlignLenFilter > 0 {
				passLen = filterByLength(i, s.AlignLenFilter)
				if !passLen {
					continue
				}
			}
			if s.FilterCigar != "" {
				passCigar = filterByCigar(i, s.FilterCigar)
				if !passCigar {
					continue
				}
			}
			if s.FilterByFlag > 0 {
				passFlag = filterByFlag(i, s.FilterByFlag)
				if !passFlag {
					continue
				}
			}
			if s.FilterByRegion != "" {
				passLocation = filterByRegion(i, bedRegion)
				if !passLocation {
					continue
				}
			}
			if s.SingleCellFormat {
				i, foundUMI = appendUmiCbc(i)
				if !foundUMI {
					continue
				}
			}
			if s.CollapseByUmi {
				newSamMap = collapseUMI(i)
				for _, j := range newSamMap {
					sam.WriteToFileHandle(out, j)
				}
			} else {
				sam.WriteToFileHandle(out, i)
			}
		}
	}
}

func usage() {
	fmt.Print("\nsamFilter -- Filters in a input SAM or BAM file based on input filters\n" +
		"Visit https://samtools.github.io/hts-specs/SAMv1.pdf for more info on SAM file specifications" +
		"Usage:\n" +
		"samFilter [options] in.sam out.sam\n")
	flag.PrintDefaults()
}

func main() {
	var alignQualityFilter *int = flag.Int("alignQualityFilter", 0, "Filter the input sam/bam file for alignments greater than the input value.")
	var alignLengthFilter *int = flag.Int("alignLengthFilter", 0, "Filter the input sam/bam for alignments longer than or equal to the input value.")
	var filterCigar *string = flag.String("filterCigar", "", "Filters the input sam/bam file for alignments that satisfy the input requirements.\nOptions:\n"+
		"\tstarrSeqIntrons - Filters out intronic alignments longer than 500 bp\n"+
		"\t\"cigarString\" - exact matches to an input CIGAR string will be kept")
	var collapseByUmi *bool = flag.Bool("collapseUMI", false, "Collapses the input sam/bam file by UMI. UMIs with highest alignment scores will kept in the event that a UMI is found multitple identical UMIs are found.\n"+
		"Must be in the gonomics single-cell annotation: \n"+
		"\tReadname appended with \"_UMI:NNNNNNNNNNNN_BX:NNNNNNNNNNNNNNNN\"\n "+
		"\tUse fastqFormat or -scFormat to append read-names with UMI and BX.")
	var scFormat *bool = flag.Bool("scFormat", false, "Appends the corrected/filtered UMI and BX to the readname from a sam/bam processed by the 10X Cellranger Count software.")
	var location *string = flag.String("coordinates", "", "Filters the input sam/bam by the input corrdinates. Only alignments that fall within the input coordinates will be kept\n"+
		"options:\n\tchr\n\tchr:start-end")
	var flagFilter *int = flag.Int("flag", 0, "Filters the input sam/bam file by SAM Flag. Only alignments with matching Flags will be kept")

	if *alignQualityFilter < 0 || *alignLengthFilter < 0 {
		log.Fatalf("the input for alingment qualtiy and length filters must be a positive intiger")
	}

	s := Settings{
		InFile:           "",
		OutFile:          "",
		AlignQualFilter:  *alignQualityFilter,
		AlignLenFilter:   *alignLengthFilter,
		FilterCigar:      *filterCigar,
		CollapseByUmi:    *collapseByUmi,
		SingleCellFormat: *scFormat,
		FilterByRegion:   *location,
		FilterByFlag:     *flagFilter,
	}

	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	runFilter(s)
}
