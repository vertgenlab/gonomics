package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"sort"
	"strings"
)

func Alt(s ScStarrSeqSettings) {
	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	var overlap []interval.Interval
	var starrSeqRead Read
	var starrSeqReadSlice []Read
	var allConstructs, umiBxSlice, allCellTypes []string
	var cluster, cbcUpdate string
	var found bool
	var multipleGems bool = false
	var gemNumber int = 1
	var ck []ClusterKey

	cellTypeMap := make(map[string]string)
	allCellTypesMap := make(map[string]int)

	if s.ScAnalysis != "" || s.CountMatrixCellTypes != "" {
		if s.ScAnalysis != "" {
			ck = ReadClusterKey(s.ScAnalysis)
		} else {
			ck = ReadClusterKey(s.CountMatrixCellTypes)
		}
		for _, i := range ck {
			allCellTypesMap[i.Cluster] = 0
			cellTypeMap[i.Bx] = i.Cluster
		}
		for i := range allCellTypesMap {
			allCellTypes = append(allCellTypes, i)
		}
		sort.Strings(allCellTypes)
	}

	bedEntries := bed.Read(s.AltMapping)
	for _, i := range bedEntries {
		tree = append(tree, i)
		allConstructs = append(allConstructs, i.Name)
	}
	selectTree = interval.BuildTree(tree)

	files := strings.Split(s.InFile, ",")
	if len(files) > 1 {
		multipleGems = true
	}
	for _, bam := range files {
		inChan, _ := sam.GoReadToChan(bam)
		for i := range inChan {
			//filter sam file for some basic quality metrics
			if !sam.FilterByQuality(i, 1) {
				continue
			}
			if !sam.FilterByCigar(i, "starrSeqIntrons") {
				continue
			}
			overlap = interval.Query(selectTree, i, "any")
			if len(overlap) != 1 {
				continue
			}
			cbc, cbcFound, _ := sam.QueryTag(i, "CB")
			umi, umiFound, _ := sam.QueryTag(i, "UB")
			if !umiFound {
				continue
			}
			if multipleGems && cbcFound {
				cbcUpdate = fmt.Sprintf("%s_%d", cbc.(string), gemNumber)
				cluster, found = cellTypeMap[cbcUpdate]
				umiBxSlice = append(umiBxSlice, fmt.Sprintf("%s_%s", umi.(string), cbcUpdate))
			} else if cbcFound {
				cbcUpdate = cbc.(string)
				cluster, found = cellTypeMap[cbcUpdate]
				umiBxSlice = append(umiBxSlice, fmt.Sprintf("%s_%s", umi.(string), cbcUpdate))
			}
			switch {
			case !cbcFound:
				starrSeqRead = Read{Bx: "NA", UMI: umi.(string), Cluster: "undefined", Construct: overlap[0].(bed.Bed).Name}
			case cbcFound && found:
				starrSeqRead = Read{Bx: cbcUpdate, UMI: umi.(string), Cluster: cluster, Construct: overlap[0].(bed.Bed).Name}
			case cbcFound && !found:
				starrSeqRead = Read{Bx: cbcUpdate, UMI: umi.(string), Cluster: "undefined", Construct: overlap[0].(bed.Bed).Name}
			}
			starrSeqReadSlice = append(starrSeqReadSlice, starrSeqRead)
		}
	}
	collapsedUmiSlice := collapseUmi(starrSeqReadSlice)
	for _, i := range collapsedUmiSlice {
		if i.Construct == "" {
			fmt.Println("collapseUmiSlice empty")
		}
	}
	r := ReadSliceAnalysisSettings{FuncSettings: s, ReadSlice: collapsedUmiSlice, AllConstructs: allConstructs, UmiBxSlice: umiBxSlice, AllCellTypes: allCellTypes,
		CellTypeMap: cellTypeMap}
	ReadSliceAnalysis(r)
}

// collapseUmi will take in a slice of Read and collapse Read entries by UMI, keeping one entry per UMI.
func collapseUmi(ssrs []Read) []Read {
	var read Read
	var found bool
	var bxUmi string
	var outSlice []Read

	umiMap := make(map[string]Read) //make a map that has cellBarcode/UMI as the key and the corresponding Read entry is the value
	for _, i := range ssrs {        // loop through slice of Read
		bxUmi = fmt.Sprintf("%s_%s", i.Bx, i.UMI) //append cell barcode and UMI together
		read, found = umiMap[bxUmi]               //check to see if the bx/umi is already in the map
		if !found {
			umiMap[bxUmi] = i //if it's not, add it with the Read to the map
		} else {
			if read.Construct == i.Construct { // if its found, check to make sure the current bx/umi maps to the same construct as the one in the map
				continue //if it does, good, do nothing. original bx/umi is retained in the map
			} else { //if the current bx/umi maps to a different construct as the construct in the map, likely a sequencing artifact somewhere. get rid of both reads.
				//could think about UMI collapsing at sam level to keep the better mapping UMI. THis is conservative, however, and doesn't seem to be throwing out too many reads.
				delete(umiMap, bxUmi)
			}
		}
	}
	for i := range umiMap {
		outSlice = append(outSlice, umiMap[i])
	}
	fmt.Printf("Found %d valid UMIs (altMapping)\n", len(outSlice))
	return outSlice
}
