package reconstruct

import (
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/bed"
	"strconv"
	"log"
)

func IlsWigToBed(postProbWig map[string]wig.Wig) [][]bed.Bed {
	keys := make([]string, len(postProbWig))

	i := 0
	for k := range postProbWig {
		keys[i] = k
		i++
	}
	// make 4 bed files
	v0Bed := make([]bed.Bed, 0)
	v1Bed := make([]bed.Bed, 0)
	v2Bed := make([]bed.Bed, 0)
	v3Bed := make([]bed.Bed, 0)
	
	ilsBeds := [][]bed.Bed{v0Bed, v1Bed, v2Bed, v3Bed}
	
	mapTopo := 0
	var maxVal float64 = 0
	if len(postProbWig["V0"].Values) != len(postProbWig["V1"].Values) || len(postProbWig["V1"].Values) != len(postProbWig["V2"].Values) || len(postProbWig["V2"].Values) != len(postProbWig["V3"].Values) {
		log.Fatal("Length of posterior probabilities not the same for each topology.")
	}
	for idx, _ := range postProbWig["V0"].Values {
		maxVal = postProbWig["V0"].Values[idx]
		mapTopo = 0
		if postProbWig["V1"].Values[idx] > maxVal {
			maxVal = postProbWig["V1"].Values[idx]
			mapTopo = 1
		} 
		
		if postProbWig["V2"].Values[idx] > maxVal{
			maxVal = postProbWig["V2"].Values[idx]
			mapTopo = 2
		} 
	
		if postProbWig["V3"].Values[idx] > maxVal{
			maxVal = postProbWig["V3"].Values[idx]
			mapTopo = 3
		}
		lastRegionIdx := len(ilsBeds[mapTopo])-1
		if len(ilsBeds[mapTopo]) == 0 { 
			ilsBeds[mapTopo] = append(ilsBeds[mapTopo], bed.Bed{Chrom: "region_0", ChromStart: idx, ChromEnd: idx, FieldsInitialized: 3})
		} else if ilsBeds[mapTopo][lastRegionIdx].GetChromEnd() == idx-1 {
			ilsBeds[mapTopo][lastRegionIdx].ChromEnd = idx
		} else {
			chromName := "region_" + strconv.Itoa(lastRegionIdx+1)
			ilsBeds[mapTopo] = append(ilsBeds[mapTopo], bed.Bed{Chrom: chromName, ChromStart: idx, ChromEnd: idx, FieldsInitialized: 3})
		}
	}

	return ilsBeds

	// pseudo logic
	// make 4 beds to return 
	// iterate through wig
		// if pos i MLE=v0/1/2/3
			// cif currently exists chrom with end position at pos i-1
				// ++chrom end position
			// else, make new chrom starting at pos i
	// return 4 beds
}