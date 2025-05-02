package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"github.com/vertgenlab/gonomics/sam"
	"golang.org/x/exp/slices"
	"log"
	"strings"
)

type BulkOutputSeqSettings struct {
	InSam     string
	InBed     string
	OutCounts string
	SingleBx  bool
	DualBx    bool
	InputNorm string
	PairedEnd bool
	RepNum    int
	ZScore    string
	CheckBx   string
}

func processPairEndSam(s BulkOutputSeqSettings, tree map[string]*interval.IntervalNode, countsMap map[string]float64) []string {
	var ansR1, ansR2 []interval.Interval
	var samBedsR1, samBedsR2 []bed.Bed
	var totPair, qualIntronFilter, twoCollapseToOne, diff, same, oneOfTwo, bothmultianddiff, none, other, totLen3_save, totLen3_discard, pairsFiltered, oneOnlyOtherFilterd, j int
	var name string

	if s.DualBx {
		countsMap = collapseDualBx(countsMap)
	}
	if s.CheckBx != "" {
		ref := fasta.ToMap(fasta.Read(s.CheckBx))
	}

	outPe := fileio.EasyCreate("/Users/sethweaver/Downloads/hsSD/starrSeq/outputseq/aln/star/testingFloor/famSS3.allreads.star.longIntron.sam")

	peReads, head := sam.GoReadSamPeToChan(s.InSam)
	sam.WriteHeaderToFileHandle(outPe, head)

	for pair := range peReads {
		totPair++

		if pair.R1.MapQ == 0 {
			qualIntronFilter++
			pair.R1 = sam.Sam{}
		} else if !filterIntrons(pair.R1) {
			qualIntronFilter++
			sam.WriteToFileHandle(outPe, pair.R1)
			pair.R1 = sam.Sam{}
		}
		if pair.R1.MapQ == 0 {
			qualIntronFilter++
			pair.R2 = sam.Sam{}
		} else if !filterIntrons(pair.R2) {
			qualIntronFilter++
			sam.WriteToFileHandle(outPe, pair.R2)
			pair.R2 = sam.Sam{}
		}
		if pair.R1.QName == "" && pair.R2.QName == "" {
			pairsFiltered++
			continue
		}

		if s.DualBx || s.SingleBx {
			ansR1 = ansR1[:0]
			ansR2 = ansR2[:0]
			if pair.R1.QName != "" {
				samBedsR1 = convert.SamToBedExcludeIntron(pair.R1)
				for j = range samBedsR1 {
					ansR1 = append(ansR1, interval.Query(tree, samBedsR1[j], "d")...)
				}
			}
			if pair.R2.QName != "" {
				samBedsR2 = convert.SamToBedExcludeIntron(pair.R2)
				for j = range samBedsR2 {
					ansR2 = append(ansR2, interval.Query(tree, samBedsR2[j], "d")...)
				}
			}
		} else {
			ansR1 = interval.Query(tree, pair.R1, "any")
			ansR2 = interval.Query(tree, pair.R2, "any")
		}

		if len(ansR1) == 2 {
			if sameConstruct(ansR1[0].(bed.Bed).Name, ansR1[1].(bed.Bed).Name, s.DualBx) {
				twoCollapseToOne++
				ansR1 = ansR1[:1]
			}
		}
		if len(ansR2) == 2 {
			if sameConstruct(ansR2[0].(bed.Bed).Name, ansR2[1].(bed.Bed).Name, s.DualBx) {
				twoCollapseToOne++
				ansR2 = ansR2[:1]
			}
		}

		switch {
		case len(ansR1) == 1 && len(ansR2) == 1:
			if sameConstruct(ansR1[0].(bed.Bed).Name, ansR2[0].(bed.Bed).Name, s.DualBx) {
				if s.DualBx {
					name = stripDualBx(ansR1[0].(bed.Bed).Name)
				} else {
					name = ansR1[0].(bed.Bed).Name
				}
				countsMap[name]++
				same++
			} else {
				//sam.WriteSamPeToFileHandle(outPe, pair)
				diff++
			}
		case len(ansR1) == 1 && len(ansR2) == 0:
			if s.DualBx {
				name = stripDualBx(ansR1[0].(bed.Bed).Name)
			} else {
				name = ansR1[0].(bed.Bed).Name
			}
			if pair.R2.QName == "" {
				oneOnlyOtherFilterd++
			}
			countsMap[name]++
			oneOfTwo++
		case len(ansR1) == 0 && len(ansR2) == 1:
			if s.DualBx {
				name = stripDualBx(ansR2[0].(bed.Bed).Name)
			} else {
				name = ansR2[0].(bed.Bed).Name
			}
			if pair.R1.QName == "" {
				oneOnlyOtherFilterd++
			}
			countsMap[name]++
			oneOfTwo++
		case (len(ansR1) == 1 && len(ansR2) == 2) || (len(ansR1) == 2 && len(ansR2) == 1):
			name = rescueTotalThree(ansR1, ansR2, s)
			if name != "" {
				countsMap[name]++
				totLen3_save++
			} else {
				totLen3_discard++
			}
		case len(ansR1) > 1 && len(ansR2) > 1:
			//fmt.Println(ansR1, ansR2)
			bothmultianddiff++
		case len(ansR1)+len(ansR2) == 0:
			none++
		default:
			//sam.WriteSamPeToFileHandle(outPe, pair)
			other++
		}
	}
	fmt.Println("Processing ", s.InSam)
	//totPair, qualIntronFilter, twoCollapseToOne, diff, same, oneOfTwo, bothmultianddiff, none, other
	fmt.Println("total read pairs: ", totPair)
	fmt.Println("reads (not pairs) filtered out by quality or intron filters: ", qualIntronFilter) //filtered out, feel good about it?
	fmt.Println("read pairs completely filtered out due to quality or intron filters: ", pairsFiltered)
	fmt.Println("reads (not pairs) that overlapped both indexes from the same construct: ", twoCollapseToOne) //counted, feel good about it
	fmt.Println("pairs in which each read maps to the same construct (opposite bx): ", same)
	fmt.Println("Read pairs in which both reads mapped to a single barcode, but to different constructs: ", diff)               // not using currently, want to see if I can rescue
	fmt.Println("Read pairs in which one read mapped to a barcode, and it's read pair didn't overlap any barcodes: ", oneOfTwo) // using currently, need to figure out why the other read doesn't map
	fmt.Println("How many of the above pairs only mapped to one total construct because the other read pair was filtered out: ", oneOnlyOtherFilterd)
	fmt.Println("Total length of three saved and discarded: ", totLen3_save, totLen3_discard)
	fmt.Println("Read pairs where both reads map to 2+ barcodes: ", bothmultianddiff) // not using need to look into
	fmt.Println("Read pairs where neither pair overlapped a barcode: ", none)         // obvi not using, haven't looked at
	fmt.Println("other, didn't match any of the other categories: ", other)

	if s.InputNorm != "" {
		inputNormalize(s.InputNorm, countsMap)
	}
	zMap := make(map[string]float64)
	if s.ZScore != "" {
		zscoreNorm(s, countsMap, zMap)
	}
	exception.PanicOnErr(outPe.Close())
	return formatResults(countsMap, zMap, s.RepNum)
}

func zscoreNorm(s BulkOutputSeqSettings, countsMap map[string]float64, zMap map[string]float64) {
	var ncCounts []float64
	nc := fileio.Read(s.ZScore)
	for i := range nc {
		ncCounts = append(ncCounts, countsMap[nc[i]])
	}
	ncAvg := numbers.AverageFloat64(ncCounts)
	ncSD := numbers.StandardDeviationFloat64(ncCounts)
	for i := range countsMap {
		zMap[i] = numbers.ZScore(countsMap[i], ncAvg, ncSD)
	}

}

func rescueTotalThree(ansR1, ansR2 []interval.Interval, s BulkOutputSeqSettings) string {
	var j int
	for i := range ansR1 {
		for j = range ansR2 {
			if sameConstruct(ansR1[i].(bed.Bed).Name, ansR2[j].(bed.Bed).Name, true) {
				if s.DualBx {
					return stripDualBx(ansR1[i].(bed.Bed).Name)
				} else {
					return ansR1[i].(bed.Bed).Name
				}
			}
		}
	}
	return ""
}

func stripDualBx(s string) string {
	slc := strings.Split(s, "_")
	return strings.Join(slc[:len(slc)-1], "_")
}

// sameConstruct returns true if the barcodes come from the same construct
func sameConstruct(c1, c2 string, dualBx bool) bool {
	if !dualBx {
		return c1 == c2
	}
	name1 := stripDualBx(c1)
	name2 := stripDualBx(c2)
	return name1 == name2
}

func formatResults(countsMap map[string]float64, zMap map[string]float64, repNum int) []string {
	var res []string

	if len(zMap) == 0 {
		for i := range countsMap {
			res = append(res, fmt.Sprintf("%s\t%.2f\t%d", i, countsMap[i], repNum))
		}
	} else {
		for i := range countsMap {
			res = append(res, fmt.Sprintf("%s\t%.2f\t%.2f\t%d", i, countsMap[i], zMap[i], repNum))
		}
	}

	slices.Sort(res)

	return res
}

func inputNormalize(normFactor string, countsMap map[string]float64) {
	var found bool
	var counts float64
	normMap := readNormFactorFile(normFactor)
	for i := range normMap {
		counts, found = countsMap[i]
		if !found {
			countsMap[i] = 0
			fmt.Printf("WARNING: %s is in the input normalization file, but no counts were found. This construct will have 0 counts in the output file\n", i)
			continue
		}
		countsMap[i] = counts * normMap[i]
	}
	checkIfEverythingWasNormed(countsMap, normMap)
}

func checkIfEverythingWasNormed(countsMap, normMap map[string]float64) {
	var found bool
	var notFound []string
	for i := range countsMap {
		_, found = normMap[i]
		if !found {
			notFound = append(notFound, i)
		}
	}
	if len(notFound) != 0 {
		log.Fatalf("ERROR: These constructs have counts but were not found in the input normalization file: %v\n", notFound)
	}
}

func makeBedTreeAndMap(inBed string) (tree map[string]*interval.IntervalNode, countsMap map[string]float64) {
	var ntrvls []interval.Interval

	countsMap = make(map[string]float64)

	bd := bed.Read(inBed)

	for i := range bd {
		ntrvls = append(ntrvls, bd[i])
		countsMap[bd[i].Name] = 0
	}

	tree = interval.BuildTree(ntrvls)

	return tree, countsMap
}

// filterIntrons returns true if they pass the intron filter
func filterIntrons(s sam.Sam) bool {
	for i := range s.Cigar {
		if s.Cigar[i].Op == 'N' && s.Cigar[i].RunLength > 700 {
			return false
		}
	}
	return true
}

func readNormFactorFile(inNorm string) map[string]float64 {
	var cols []string
	normFactorMap := make(map[string]float64)
	in := fileio.Read(inNorm)
	for i := 1; i < len(in); i++ {
		cols = strings.Split(in[i], "\t")
		normFactorMap[cols[0]] = parse.StringToFloat64(cols[3])
	}
	return normFactorMap
}

func checkBx(s sam.Sam, bxCoord bed.Bed, ref map[string][]fasta.Fasta) {
	sam.SamBedToBases(s, bxCoord)
}

// non-pe
func BulkOutputSeq(s BulkOutputSeqSettings) []string {
	var ans []interval.Interval
	var zero, multi int

	tree, countsMap := makeBedTreeAndMap(s.InBed)

	if s.PairedEnd {
		return processPairEndSam(s, tree, countsMap)
	}

	inSam, _ := sam.GoReadToChan(s.InSam)

	for read := range inSam {
		if read.MapQ < 1 {
			continue
		} else if !filterIntrons(read) {
			continue
		}
		if s.DualBx || s.SingleBx {
			ans = interval.Query(tree, read, "d")
		} else {
			ans = interval.Query(tree, read, "any")
		}

		if len(ans) == 2 {
			if sameConstruct(ans[0].(bed.Bed).Name, ans[1].(bed.Bed).Name, s.DualBx) {
				ans = ans[:1]
			}
		}

		switch len(ans) {
		case 0:
			zero++
			continue
		case 1:
			countsMap[ans[0].(bed.Bed).Name]++
		default:
			multi++
			continue
		}
	}
	if s.DualBx {
		countsMap = collapseDualBx(countsMap)
	}
	if s.InputNorm != "" {
		inputNormalize(s.InputNorm, countsMap)
	}
	zMap := make(map[string]float64)
	if s.ZScore != "" {
		zscoreNorm(s, countsMap, zMap)
	}
	return formatResults(countsMap, zMap, s.RepNum)
}
