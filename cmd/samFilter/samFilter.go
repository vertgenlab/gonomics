package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
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
	SortByPosition   bool
	OutBam           bool
	NoHeader         bool
}

func filterByQuality(a sam.Sam, filter int) bool {
	if int(a.MapQ) >= filter {
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

func collapseUMI(inMap map[string]sam.Sam, inRead sam.Sam) map[string]sam.Sam {
	umiMap := inMap
	sc := sam.ToSingleCellAlignment(inRead)
	umiBxString := fmt.Sprintf("%s%s", dna.BasesToString(sc.Bx), dna.BasesToString(sc.Umi))
	_, found := umiMap[umiBxString]
	if !found { //if the UMI isn't in the map, add it.
		umiMap[umiBxString] = inRead
	} else {
		if umiMap[umiBxString].MapQ > inRead.MapQ {
			//if the UMI is in the map and the UMI in the map has a better alignment score -- do nothing.
		} else if umiMap[umiBxString].MapQ == inRead.MapQ {
			if umiMap[umiBxString].Pos >= inRead.Pos { //if the UMIs are tied on alignment score, prioritize the one with a higher start position since GWAS STARR-seq construct barcodes are on the end of constructs,
				//higher chance that the read contains the construct barcode.
				//do nothing
			} else if umiMap[umiBxString].Pos < inRead.Pos { //if the UMI in the map has a lower position than the query, replace it.
				umiMap[umiBxString] = inRead
			}
		} else if umiMap[umiBxString].MapQ < inRead.MapQ { //if the UMI in the map has a lower alignment score than the query, replace it.
			umiMap[umiBxString] = inRead
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
	var umiMap = make(map[string]sam.Sam)
	var words, p []string
	var chrom string
	var start, end int
	var bedRegion bed.Bed
	var outSlice []sam.Sam
	var err error
	var bw *sam.BamWriter

	inChan, header := sam.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)
	if s.OutBam {
		bw = sam.NewBamWriter(out, header)
	} else if s.OutFile != "stdout" && s.NoHeader == false {
		sam.WriteHeaderToFileHandle(out, header)
	}
	if s.FilterByRegion != "" {
		words = strings.Split(s.FilterByRegion, ":")
		chrom = words[0]
		if len(words) == 1 {
			bedRegion = bed.Bed{Chrom: chrom, ChromStart: 0, ChromEnd: numbers.MaxInt}
		} else {
			p = strings.Split(words[1], "-")
			start = common.StringToInt(p[0])
			end = common.StringToInt(p[1])
			bedRegion = bed.Bed{Chrom: chrom, ChromStart: start, ChromEnd: end}
		}
	}
	for i := range inChan {
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
			umiMap = collapseUMI(umiMap, i)
			continue
		} else if s.SortByPosition {
			outSlice = append(outSlice, i)
		} else {
			if s.OutBam {
				sam.WriteToBamFileHandle(bw, i, 0)
			} else {
				sam.WriteToFileHandle(out, i)
			}
		}
	}
	if s.CollapseByUmi {
		if s.SortByPosition {
			for _, i := range umiMap {
				outSlice = append(outSlice, i)
			}
		} else {
			for _, i := range umiMap {
				if s.OutBam {
					sam.WriteToBamFileHandle(bw, i, 0)
				} else {
					sam.WriteToFileHandle(out, i)
				}
			}
		}
	}
	if s.SortByPosition {
		sam.SortByCoord(outSlice)
		for _, i := range outSlice {
			if s.OutBam {
				sam.WriteToBamFileHandle(bw, i, 0)
			} else {
				sam.WriteToFileHandle(out, i)
			}
		}
	}
	if s.OutBam {
		err = bw.Close()
		exception.PanicOnErr(err)
	} else {
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print("\nsamFilter -- Filters in a input SAM or BAM file based on input filters\n" +
		"Visit https://samtools.github.io/hts-specs/SAMv1.pdf for more info on SAM file specifications\n" +
		"The output filename can be replaced with stdout and the sam file without a header will be printed to the screen.\n\n" +
		"Usage:\n" +
		"samFilter [options] in.sam out.sam\n\noptions:\n")
	flag.PrintDefaults()
}

func main() {
	var alignQualityFilter *int = flag.Int("alignQualityFilter", 0, "Filter the input sam/bam file for alignments with MAPQ equal or greater than the input value.")
	var alignLengthFilter *int = flag.Int("alignLengthFilter", 0, "Filter the input sam/bam for alignments longer than or equal to the input value.")
	var filterCigar *string = flag.String("filterCigar", "", "Filters the input sam/bam file for alignments that satisfy the input requirements.\nOptions:\n"+
		"\tstarrSeqIntrons - Filters out intronic alignments longer than 500 bp\n"+
		"\t\"cigarString\" - exact matches to an input CIGAR string will be kept")
	var collapseByUmi *bool = flag.Bool("collapseUMI", false, "Collapses the input sam/bam file by UMI. UMIs with highest alignment scores will kept in the event that a UMI is found multiple identical UMIs are found.\n"+
		"Must be in the gonomics single-cell annotation: \n"+
		"\tReadname appended with \"_UMI:NNNNNNNNNNNN_BX:NNNNNNNNNNNNNNNN\"\n "+
		"\tUse fastqFormat or -scFormat to append read-names with UMI and BX.")
	var scFormat *bool = flag.Bool("scFormat", false, "Appends the corrected/filtered UMI and BX to the readname from a sam/bam processed by the 10X Cellranger Count software. "+
		"Cell barcode -- CB tag. UMI -- UB tag.  ***NOTE: scFormat is only compatible with BAM input")
	var location *string = flag.String("coordinates", "", "Filters the input sam/bam by the input coordinates. Only alignments that fall within the input coordinates will be kept\n"+
		"Format options:\n\tchr\n\tchr:start-end")
	var flagFilter *int = flag.Int("flag", 0, "Filters the input sam/bam file by SAM Flag. Only alignments with matching Flags will be kept")
	var sortByPosition *bool = flag.Bool("sort", false, "Sorts the output sam/bam file by position")
	var outBam *bool = flag.Bool("bam", false, "Write the output in BAM format")
	var noHeader *bool = flag.Bool("noHeader", false, "Sam header isn't included in the output file. This option is only compatible with the Sam output type")

	if *alignQualityFilter < 0 || *alignLengthFilter < 0 {
		log.Fatalf("the input for alignment qualtiy and length filters must be a positive intiger")
	}

	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	if strings.HasSuffix(flag.Arg(0), ".sam") && *scFormat {
		log.Fatalf("-scFormat is only compatable with a BAM input")
	}
	s := Settings{
		InFile:           flag.Arg(0),
		OutFile:          flag.Arg(1),
		AlignQualFilter:  *alignQualityFilter,
		AlignLenFilter:   *alignLengthFilter,
		FilterCigar:      *filterCigar,
		CollapseByUmi:    *collapseByUmi,
		SingleCellFormat: *scFormat,
		FilterByRegion:   *location,
		FilterByFlag:     *flagFilter,
		SortByPosition:   *sortByPosition,
		OutBam:           *outBam,
		NoHeader:         *noHeader,
	}
	runFilter(s)
}
