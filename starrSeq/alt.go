package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
)

func Alt(s ScStarrSeqSettings) {
	var tree = make([]interval.Interval, 0)
	var selectTree map[string]*interval.IntervalNode
	var overlap []interval.Interval
	var starrSeqRead Read
	var starrSeqReadSlice []Read

	out := fileio.EasyCreate(s.OutFile)

	bedEntries := bed.Read(s.AltMapping)
	for _, i := range bedEntries {
		tree = append(tree, i)
	}
	selectTree = interval.BuildTree(tree)

	inChan, _ := sam.GoReadToChan(s.InFile)

	for i := range inChan {
		overlap = interval.Query(selectTree, i, "any")
		if len(overlap) != 1 {
			continue
		}
		cbc, _, _ := sam.QueryTag(i, "CB")
		umi, _, _ := sam.QueryTag(i, "UB")
		starrSeqRead = Read{Bx: cbc.(string), UMI: umi.(string), Construct: overlap[0].(bed.Bed).Name}
		starrSeqReadSlice = append(starrSeqReadSlice, starrSeqRead)
	}
	collapsedUmiSlice := collapseUmi(starrSeqReadSlice)
	constructMap := readSliceToPseudobulk(s, collapsedUmiSlice)
	WritePseudobulkMap(constructMap, out)
	err := out.Close()
	exception.PanicOnErr(err)
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

func readSliceToPseudobulk(s ScStarrSeqSettings, readSlice []Read) map[string]float64 {
	constructMap := make(map[string]float64)
	beds := bed.Read(s.AltMapping)
	for _, i := range beds {
		constructMap[i.Name] = 0
	}
	for _, i := range readSlice {
		constructMap[i.Construct] += 1
	}
	if s.InputNormalize != "" {
		InputNormalize(constructMap, s.InputNormalize)
	}
	return constructMap
}
