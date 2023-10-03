package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"sort"
)

// ParseInputSequencingSam reads and parses an alignment file for an input library sequencing run. It creates a read-count map for all constructs in the provided GTF
func ParseInputSequencingSam(s ScStarrSeqSettings) {
	var tree = make([]interval.Interval, 0)
	var bedSizeMap = make(map[string]int)
	var countsMap = make(map[string]float64)
	var selectTree map[string]*interval.IntervalNode
	var constructName string
	var counts float64
	var currEntry sam.Sam

	inChan, _ := sam.GoReadToChan(s.InFile)
	bedEntries := bed.Read(s.InputSequencing)

	for _, i := range bedEntries {
		tree = append(tree, i)
		bedSizeMap[i.Name] = bed.Size(i)
		countsMap[i.Name] = 0
	}
	selectTree = interval.BuildTree(tree)

	currEntry = <-inChan
	for i := range inChan {
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
	//normalize reads to 500bp
	for i := range countsMap {
		counts, _ = countsMap[i]
		countsMap[i] = counts / (float64(bedSizeMap[i]) / 500.0)
	}
	calculateNormFactor(countsMap, s.OutFile)
}

// calculateNormFactor takes a countsMap created in ReadInputSequencingSam and writes out a data frame with the name of the construct, counts per 500bp, percent abundance, and normalization factor
func calculateNormFactor(countsMap map[string]float64, outFile string) {
	var totalReads float64 = 0
	var entry string
	var percentLib float64
	var matrix []string

	out := fileio.EasyCreate(outFile)
	fileio.WriteToFileHandle(out, "construct\treadsPer500bp\tpercentLibrary\tnormFactor")

	numConstructs := len(countsMap)
	idealPerc := 100.0 / float64(numConstructs)

	for i := range countsMap {
		totalReads = countsMap[i] + totalReads
	}

	for i := range countsMap {
		percentLib = (countsMap[i] / totalReads) * 100
		entry = fmt.Sprintf("%s\t%f\t%f\t%f", i, countsMap[i], percentLib, idealPerc/percentLib)
		matrix = append(matrix, entry)
	}
	sort.Strings(matrix)
	for _, i := range matrix {
		fileio.WriteToFileHandle(out, i)
	}
	err := out.Close()
	exception.PanicOnErr(err)
}
