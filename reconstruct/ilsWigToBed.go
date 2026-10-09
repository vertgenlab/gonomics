package reconstruct

import (
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/bed"
	"log"
)

func IlsWigToBed(postProbWig map[string]wig.Wig) [][]bed.Bed {
	keys := make([]string, len(postProbWig))

	i := 0
	for k := range postProbWig {
		keys[i] = k
		i++
	}

	// make 4 bed files and output
	ilsBeds := make([][]bed.Bed, 4)
	trackNames := []string{"V0", "V1", "V2", "V3"}

	n := len(postProbWig[trackNames[0]].Values)
    for _, name := range trackNames[1:] {
        if len(postProbWig[name].Values) != n {
            log.Fatal("Length of posterior probabilities not the same for each topology.")
        }
    }

	w0 := postProbWig[trackNames[0]]
	start := w0.Start
	step := w0.Step
	chromName := "chr1"
	var maxVal float64 = 0
	var maxTopo int = 0

	val := postProbWig[trackNames[0]].Values[0]
	var bedIdx int = 0
	for idx := 0; idx < n; idx++ {
		maxVal = postProbWig[trackNames[0]].Values[idx]
		maxTopo = 0
		for t := 1; t < 4; t++ {
			val = postProbWig[trackNames[t]].Values[idx]
			if val > maxVal {
				maxVal = val
				maxTopo = t
			}
		}

		lastRegionIdx := len(ilsBeds[maxTopo])-1
		bedIdx = (start + idx*step) - 1
		if len(ilsBeds[maxTopo]) == 0 { //// start new bed region run, no existing regions in bed
			ilsBeds[maxTopo] = append(ilsBeds[maxTopo], bed.Bed{Chrom: chromName, ChromStart: bedIdx, ChromEnd: bedIdx+1, Name: trackNames[maxTopo], FieldsInitialized: 4})
		} else {
			if ilsBeds[maxTopo][lastRegionIdx].GetChromEnd() == bedIdx { //// extend existing bed region by 1
			ilsBeds[maxTopo][lastRegionIdx].ChromEnd = bedIdx+1
			} else { //// start new bed region run
				ilsBeds[maxTopo] = append(ilsBeds[maxTopo], bed.Bed{Chrom: chromName, ChromStart: bedIdx, ChromEnd: bedIdx+1, Name: trackNames[maxTopo], FieldsInitialized: 4})
			}
		}
	}

	return ilsBeds

}