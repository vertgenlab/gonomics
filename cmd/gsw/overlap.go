package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/fileio"
	"sync"
	"strings"

)

func overlapSelect(selectTarget string, bedFmt string, output string, invert bool) {
	switch true {
	case strings.HasSuffix(selectTarget, ".axt"):
		findAxtBedOverlap(selectTarget, bedFmt, output, invert)

	}
}

func findChainBedOverlap(chainFmt string, bedFmt string, output string, invert bool) {

}

//axt are select/target regions we are looking for overlaps from
func findAxtBedOverlap(axtFmt string, bedFmt string, output string, invert bool) {
	target := mkAxtMap(axtFmt)
	query := make(chan *bed.Bed)
	
	go bed.ReadToChan(bedFmt, query)

	var wg sync.WaitGroup
	wg.Add(1)

	go overlapAxtBed(output, target, query, invert, &wg)
	wg.Wait()
}

func overlapAxtBed(output string, target map[string][]*axt.Axt, query <-chan *bed.Bed, invert bool, wg *sync.WaitGroup) {
	answer := fileio.MustCreate(output)
	defer answer.Close()
	for eachRegion := range query {
		if overlapAxtCompareBeds(target[eachRegion.Chrom], eachRegion) {
			//TODO: consider changing this so people can input n fields
			bed.WriteBed(answer, eachRegion, 5)
		} else {
			if invert {
				bed.WriteBed(answer, eachRegion, 5)
			}
		}
	}
	wg.Done()
}

//check if a single bed record has overlap
func overlapAxtCompareBeds(target []*axt.Axt, query *bed.Bed) bool {
	for _, each := range target {
		if axt.OverlapAxtBed(each, query) {
			return true
		}
	}
	return false
}

func mkAxtMap(axtFmt string) map[string][]*axt.Axt {
	input := fileio.EasyOpen(axtFmt)
	reader := make(chan *axt.Axt)

	target := make(map[string][]*axt.Axt)

	go axt.ReadToChan(input, reader)
	for each := range reader {
		target[each.RName] = append(target[each.RName], each)
	}
	return target
}

//chain bed overlap
func mkChainMap(chainFmt string) map[string][]*chain.Chain {
	input := fileio.EasyOpen(chainFmt)
	reader := make(chan *chain.Chain)

	target := make(map[string][]*chain.Chain)

	go chain.ReadToChan(input, reader)

	for each := range reader {
		target[each.TName] = append(target[each.TName], each)
	}
	return target
}

