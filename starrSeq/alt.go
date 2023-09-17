package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

func Alt(s ScStarrSeqSettings) {
	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	var overlap []interval.Interval
	var starrSeqRead Read
	var starrSeqReadSlice []Read
	var allConstructs, allCellTypes, umiBxSlice []string
	var cluster, cbcUpdate string
	var found bool
	var multipleGems bool = false
	var gemNumber int = 1

	cellTypeMap := make(map[string]string)

	if s.ScAnalysis != "" {
		ck := ReadClusterKey(s.ScAnalysis)
		for _, i := range ck {
			allCellTypes = append(allConstructs, i.Cluster)
			cellTypeMap[i.Bx] = i.Cluster
		}
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
			overlap = interval.Query(selectTree, i, "any")
			if len(overlap) != 1 {
				continue
			}
			cbc, tagFound, _ := sam.QueryTag(i, "CB")
			umi, _, _ := sam.QueryTag(i, "UB")
			if multipleGems && tagFound {
				cbcUpdate = fmt.Sprintf("%s_%d", cbc.(string), gemNumber)
				cluster, found = cellTypeMap[cbcUpdate]
				umiBxSlice = append(umiBxSlice, fmt.Sprintf("%s_%s", umi.(string), cbcUpdate))
			} else if tagFound {
				cbcUpdate = cbc.(string)
				cluster, found = cellTypeMap[cbcUpdate]
				umiBxSlice = append(umiBxSlice, fmt.Sprintf("%s_%s", umi.(string), cbcUpdate))
			}
			switch {
			case !tagFound:
				starrSeqRead = Read{Bx: "NA", UMI: umi.(string), Cluster: "undefined", Construct: overlap[0].(bed.Bed).Name}
			case tagFound && found:
				starrSeqRead = Read{Bx: cbcUpdate, UMI: umi.(string), Cluster: cluster, Construct: overlap[0].(bed.Bed).Name}
			case tagFound && !found:
				starrSeqRead = Read{Bx: cbcUpdate, UMI: umi.(string), Cluster: "undefined", Construct: overlap[0].(bed.Bed).Name}
			}
			starrSeqReadSlice = append(starrSeqReadSlice, starrSeqRead)
		}
	}
	collapsedUmiSlice := collapseUmi(starrSeqReadSlice)
	r := ReadSliceAnalysisSettings{FuncSettings: s, ReadSlice: collapsedUmiSlice, AllConstructs: allConstructs, UmiBxSlice: umiBxSlice, AllCellTypes: allCellTypes,
		CellTypeMap: cellTypeMap}
	ReadSliceAnalysis(r)
}

func collapseUmi(ssrs []Read) []Read {
	var read Read
	var found bool
	var bxUmi string
	var outSlice []Read

	umiMap := make(map[string]Read)
	for _, i := range ssrs {
		bxUmi = fmt.Sprintf("%s_%s", i.Bx, i.UMI)
		read, found = umiMap[bxUmi]
		if !found {
			umiMap[bxUmi] = i
		} else {
			if read.Construct == i.Construct {
				continue
			} else {
				umiMap[bxUmi] = Read{
					Bx:        "",
					Construct: "",
					UMI:       "",
					Cluster:   "",
				}
			}
		}
	}
	for i := range umiMap {
		outSlice = append(outSlice, umiMap[i])
	}
	fmt.Printf("Found %d valid UMIs (altMapping)\n", len(outSlice))
	return outSlice
}
