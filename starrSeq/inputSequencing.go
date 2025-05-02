package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"sort"
	"strings"
)

type InputSeqSettings struct {
	InSam         string
	InBed         string
	Outfile       string
	DualBx        bool
	Plasmidsaurus bool
	PairedEnd     bool
}

// ParseInputSequencingSam reads and parses an alignment file for an input library sequencing run. It creates a read-count map for all constructs in the provided bed file.
// input sam file must be position sorted
func ParseInputSequencingSam(s InputSeqSettings) {
	var tree = make([]interval.Interval, 0)
	var bedSizeMap = make(map[string]int)
	var countsMap = make(map[string]float64)
	var selectTree map[string]*interval.IntervalNode
	var constructName string
	var counts float64
	var currEntry sam.Sam

	bedEntries := bed.Read(s.InBed)

	for _, i := range bedEntries {
		tree = append(tree, i)
		bedSizeMap[i.Name] = bed.Size(i)
		countsMap[i.Name] = 0
	}
	selectTree = interval.BuildTree(tree)

	if s.PairedEnd {
		pairedEndMode(s, selectTree, countsMap, bedSizeMap)
	}

	inChan, _ := sam.GoReadToChan(s.InSam)

	if s.Plasmidsaurus {
		plasmidsaurus(s, inChan, selectTree, countsMap)
		return
	}

	currEntry = <-inChan
	for i := range inChan {
		if !sam.FilterByQuality(i, 0) { //filter out MAPQ = 0 reads
			continue
		}
		//remove PCR duplicates
		if interval.AreEqual(currEntry, i) {
			if currEntry.Qual >= i.Qual {
				continue
			} else {
				currEntry = i
				continue
			}
		}
		overlappedBedEntry := interval.Query(selectTree, currEntry, "any")
		if len(overlappedBedEntry) != 1 {
			currEntry = i
			continue
		}
		constructName = overlappedBedEntry[0].(bed.Bed).Name
		counts, _ = countsMap[constructName]
		countsMap[constructName] = counts + 1
		currEntry = i
	}
	//do the loop for the last read in the file
	overlappedBedEntry := interval.Query(selectTree, currEntry, "any")
	if len(overlappedBedEntry) == 1 {
		constructName = overlappedBedEntry[0].(bed.Bed).Name
		counts, _ = countsMap[constructName]
		countsMap[constructName] = counts + 1
	}

	if s.DualBx {
		countsMap = collapseDualBx(countsMap)
	}

	if !s.DualBx {
		//normalize reads to 500bp
		for i := range countsMap {
			counts, _ = countsMap[i]
			countsMap[i] = counts / (float64(bedSizeMap[i]) / 500.0)
		}
	}
	calculateNormFactor(s, countsMap)
}

func pairedEndMode(s InputSeqSettings, bedTree map[string]*interval.IntervalNode, countsMap map[string]float64, bedSizeMap map[string]int) {
	var ansR1, ansR2 []interval.Interval
	var none, ideal, oneOverlapEachDiff, tot, passFilter, onlyOne, onlyOneBcFilter, other int
	var counts float64

	if s.DualBx {
		countsMap = collapseDualBx(countsMap)
	}

	readsChan, _ := sam.GoReadSamPeToChan(s.InSam)

	for pair := range readsChan {
		tot++
		if pair.R1.MapQ == 0 {
			pair.R1 = sam.Sam{}
		}
		if pair.R2.MapQ == 0 {
			pair.R2 = sam.Sam{}
		}
		if pair.R1.QName == "" && pair.R2.QName == "" {
			continue
		}
		passFilter++
		if s.DualBx {
			ansR1 = interval.Query(bedTree, pair.R1, "d")
			ansR2 = interval.Query(bedTree, pair.R2, "d")
		} else {
			ansR1 = interval.Query(bedTree, pair.R1, "any")
			ansR2 = interval.Query(bedTree, pair.R2, "any")
		}

		switch {
		case len(ansR1)+len(ansR2) == 0:
			none++
			continue
		case len(ansR1) == 1 && len(ansR2) == 1:
			if sameConstruct(ansR1[0].(bed.Bed).Name, ansR2[0].(bed.Bed).Name, s.DualBx) {
				ideal++
				if s.DualBx {
					countsMap[stripDualBx(ansR1[0].(bed.Bed).Name)]++
				} else {
					countsMap[ansR1[0].(bed.Bed).Name]++
				}
			} else {
				oneOverlapEachDiff++
			}
		case len(ansR1)+len(ansR2) == 1:
			onlyOne++
			if pair.R1.QName == "" || pair.R2.QName == "" {
				onlyOneBcFilter++
			}
			if len(ansR1) == 1 {
				if s.DualBx {
					countsMap[stripDualBx(ansR1[0].(bed.Bed).Name)]++
				} else {
					countsMap[ansR1[0].(bed.Bed).Name]++
				}
			} else {
				if s.DualBx {
					countsMap[stripDualBx(ansR2[0].(bed.Bed).Name)]++
				} else {
					countsMap[ansR2[0].(bed.Bed).Name]++
				}
			}
		default:
			other++
		}

	}
	fmt.Println("Total pairs analysed: ", tot)
	fmt.Println("Pairs passing filters: ", passFilter)
	fmt.Println("No overlaps at all: ", none)
	fmt.Println("Ideal case: ", ideal)
	fmt.Println("Each one overlap, but diff constructs: ", oneOverlapEachDiff)
	fmt.Println("one total overlap: ", onlyOne)
	fmt.Println("Only one total overlap because the other read in the pair was filtered: ", onlyOneBcFilter)
	fmt.Println("other cases: ", other)

	if !s.DualBx {
		//normalize reads to 500bp
		for i := range countsMap {
			counts, _ = countsMap[i]
			countsMap[i] = counts / (float64(bedSizeMap[i]) / 500.0)
		}
	}
	calculateNormFactor(s, countsMap)

}

func collapseDualBx(countsMap map[string]float64) map[string]float64 {
	var name string
	var slc []string
	collapseMap := make(map[string]float64)

	for i := range countsMap {
		slc = strings.Split(i, "_")
		name = strings.Join(slc[:len(slc)-1], "_")
		collapseMap[name] += countsMap[i]
	}
	if len(collapseMap)*2 != len(countsMap) {
		fmt.Printf("WARNING: Dual barcode collapsed map does not have twice as fewer constructs as the original map, some constructs may not have collapsed.")
	}
	return collapseMap
}

// calculateNormFactor takes a countsMap created in ReadInputSequencingSam and writes out a data frame with the name of the construct, counts per 500bp, percent abundance, and normalization factor
func calculateNormFactor(s InputSeqSettings, countsMap map[string]float64) {
	var totalReads float64 = 0
	var entry string
	var percentLib float64
	var matrix []string

	out := fileio.EasyCreate(s.Outfile)
	if s.DualBx {
		fileio.WriteToFileHandle(out, "construct\treadCounts\tpercentLibrary\tnormFactor")
	} else {
		fileio.WriteToFileHandle(out, "construct\treadsPer500bp\tpercentLibrary\tnormFactor")
	}

	numConstructs := len(countsMap)
	idealPerc := 100.0 / float64(numConstructs)

	for i := range countsMap {
		totalReads = countsMap[i] + totalReads
	}

	for i := range countsMap {
		percentLib = (countsMap[i] / totalReads) * 100
		entry = fmt.Sprintf("%s\t%.0f\t%f\t%f", i, countsMap[i], percentLib, idealPerc/percentLib)
		matrix = append(matrix, entry)
	}
	sort.Strings(matrix)
	for _, i := range matrix {
		fileio.WriteToFileHandle(out, i)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}

/*
func plamsidsaurus(s InputSeqSettings, samChan <-chan sam.Sam, bedTree map[string]*interval.IntervalNode, countsMap map[string]float64) {
	var ans []interval.Interval
	var none, multiple int
	var found bool
	bxMap := make(map[string][]bed.Bed)
	bx := bed.Read("/Users/sethweaver/Downloads/cdkn1aMUTATIONstarr/inputSeq/ref/cdkn1a_mut.bx.bed")

	for i := range bx {
		if !strings.Contains(bx[i].Name, "CDKN1A") {
			continue
		}
		_, found = bxMap[bx[i].Name[:len(bx[i].Name)-2]]
		if !found {
			bxMap[bx[i].Name[:len(bx[i].Name)-2]] = []bed.Bed{bx[i]}
		} else {
			bxMap[bx[i].Name[:len(bx[i].Name)-2]] = append(bxMap[bx[i].Name[:len(bx[i].Name)-2]], bx[i])
		}
	}

	var a, b, c, d, oor int

	for i := range samChan {
		ans = interval.Query(bedTree, i, "any")
		switch len(ans) {
		case 0:
			none++
		case 1:
			_, found = bxMap[ans[0].(bed.Bed).Name]
			if found {
				a, b, c, d, oor = doWork(i, bxMap[ans[0].(bed.Bed).Name], a, b, c, d, oor)
			}
		default:
			multiple++

		}
	}
	fmt.Println(a, b, c, d, oor)
	fmt.Printf("None: %d\nMulitple: %d\n", none, multiple)
}

func doWork(s sam.Sam, bx []bed.Bed, a, b, c, d, oor int) (int, int, int, int, int) {
	var bases []dna.Base
	var start, end int
	if interval.Within(bx[0], s) {
		bases = dna.StringToBases(strings.Split(bx[0].Name, "_")[len(strings.Split(bx[0].Name, "_"))-2])
		start, end = bedToSamIdx(s, bx[0])
		if end > len(s.Seq) {
			oor++
		} else {
			if dna.CompareSeqsIgnoreCase(bases, s.Seq[start:end]) == 0 {
				a++
			} else {
				fmt.Println(dna.BasesToString(bases), dna.BasesToString(s.Seq[start:end]))
				b++
			}
		}
	}
	if interval.Within(bx[1], s) {
		bases = dna.StringToBases(strings.Split(bx[1].Name, "_")[len(strings.Split(bx[0].Name, "_"))-2])
		start, end = bedToSamIdx(s, bx[1])
		if dna.CompareSeqsIgnoreCase(bases, s.Seq[start:end]) == 0 {
			c++
		} else {
			d++
		}
	}
	return a, b, c, d, oor
}

func bedToSamIdx(s sam.Sam, b bed.Bed) (int, int) {
	var start, end int
	start = b.ChromStart - s.GetChromStart()
	end = start + bed.Size(b)
	return start, end
}

*/

func plasmidsaurus(s InputSeqSettings, samChan <-chan sam.Sam, bedTree map[string]*interval.IntervalNode, countsMap map[string]float64) {
	var ans []interval.Interval
	var none, multiple, supp, good int

	out := fileio.EasyCreate(s.Outfile)

	for i := range samChan {
		ans = interval.Query(bedTree, i, "any")
		switch len(ans) {
		case 0:
			none++
		case 1:
			if i.Flag&2048 == 2048 || i.MapQ == 0 {
				supp++
			} else {
				countsMap[ans[0].(bed.Bed).Name]++
				good++
			}
		default:
			multiple++
		}
	}
	for i := range countsMap {
		fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%f", i, countsMap[i]))
	}
	exception.PanicOnErr(out.Close())
}
