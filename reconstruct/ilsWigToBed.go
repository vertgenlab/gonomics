package reconstruct

import (
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/wig"
)

// IlsWigToBed takes in n posterior probability wigs and converts to n beds, one for each topology
// by taking the map estimate over these posterior probabilities for each position
// and assigning that position to the correct topology's bed
// e.g. for n=4, returns a slice of 4 beds containing info on which positions are in
// which topology states V0-V3
// chromName defaults to chr1
func IlsWigToBed(postProbWig map[string]wig.Wig, chromName string) [][]bed.Bed {
	// check number of topologies
	keys := make([]string, len(postProbWig))

	numTopos := 0
	for topoIdx := range postProbWig {
		keys[numTopos] = topoIdx
		numTopos++
	}

	if numTopos == 0 {
		log.Fatal("Must have at least 1 wig track.")
	}

	// make 4 bed files and output template
	ilsBeds := make([][]bed.Bed, numTopos)
	trackNames := make([]string, numTopos)
	for topoIdx := range trackNames {
		trackNames[topoIdx] = fmt.Sprintf("V%d", topoIdx)
	}

	// check all wig tracks cover same regions
	n := len(postProbWig[trackNames[0]].Values)
	trackStart := postProbWig[trackNames[0]].Start
	for _, name := range trackNames[1:] {
		if len(postProbWig[name].Values) != n {
			log.Fatal("Length of posterior probabilities not the same for each topology.")
		}

		if postProbWig[name].Start != trackStart {
			log.Fatal("Tracks do not start at the same position.")
		}

		if postProbWig[name].Step != 1 {
			log.Fatal("Tracks must have a step of 1.")
		}
	}

	// initialise
	w0 := postProbWig[trackNames[0]]
	start := w0.Start
	step := w0.Step
	if chromName == "" {
		chromName = "chr1"
	}
	var maxVal float64 = 0
	var maxTopo int = 0

	val := postProbWig[trackNames[0]].Values[0]
	var bedIdx int = 0
	for idx := 0; idx < n; idx++ {
		maxVal = postProbWig[trackNames[0]].Values[idx]
		maxTopo = 0
		for t := 1; t < numTopos; t++ {
			val = postProbWig[trackNames[t]].Values[idx]
			if val > maxVal {
				maxVal = val
				maxTopo = t
			}
		}

		lastRegionIdx := len(ilsBeds[maxTopo]) - 1
		bedIdx = (start + idx*step) - 1
		if len(ilsBeds[maxTopo]) == 0 { //// start new bed region run, no existing regions in bed
			ilsBeds[maxTopo] = append(ilsBeds[maxTopo], bed.Bed{Chrom: chromName, ChromStart: bedIdx, ChromEnd: bedIdx + 1, Name: trackNames[maxTopo], FieldsInitialized: 4})
		} else {
			if ilsBeds[maxTopo][lastRegionIdx].GetChromEnd() == bedIdx { //// extend existing bed region by 1
				ilsBeds[maxTopo][lastRegionIdx].ChromEnd = bedIdx + 1
			} else { //// start new bed region run
				ilsBeds[maxTopo] = append(ilsBeds[maxTopo], bed.Bed{Chrom: chromName, ChromStart: bedIdx, ChromEnd: bedIdx + 1, Name: trackNames[maxTopo], FieldsInitialized: 4})
			}
		}
	}

	return ilsBeds

}
