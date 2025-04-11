package reconstruct

import (
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"encoding/csv"
	"fmt"
	"log"
	"strconv"
	"os"

)

func readPostProbs(postProbsFilename string) map[string][]float64 {
	f, err := os.Open(postProbsFilename)
    if err != nil {
        log.Fatal("Unable to read input file " + postProbsFilename, err)
    }
    defer f.Close()

	reader := csv.NewReader(f)
	mapPostProbs := make(map[string][]float64)
	
	for {
		row, err := reader.Read()
		if err != nil {
			break
		}
		if len(row) != 5 {
			// hardcoded to ignore anything that isn't n_abc=3,n_ab=3 for 4 different topologies
			continue
		}

		probs := make([]float64, 4)
		for idx := 1; idx <= 4; idx++ {
			probs[idx-1], err = strconv.ParseFloat(row[idx], 64)
			if err != nil {
				log.Fatalf("Couldn't convert float %s", row[idx])
			}
		}

		mapPostProbs[row[0]] = probs
	}

	return mapPostProbs
}

func PostProbToWig(postProbsFile string, mafInput []*maf.Maf) map[string]wig.Wig {
	
	if len(mafInput) != 1 {
		log.Fatal("Assumes only 1 maf block")
	}

	prevLength := 0
	curLength := 0
	for idx, species := range mafInput[0].Species {
		curLength = species.SLine.Size
		if idx == 0 {
			prevLength = curLength
		} else if prevLength != curLength {
			log.Fatal("Species must have same sequence length.")
		}
	}
	
	chroms := make([]chromInfo.ChromInfo, 0)
	var trackName string
	for idx := 0; idx <= 3; idx++ {
		trackName = fmt.Sprintf("%s%d", "V", idx)
		chroms = append(chroms, chromInfo.ChromInfo{Name: trackName, Size: curLength, Order: idx})
	}
	refChromSizes := chromInfo.SliceToMap(chroms)
	outputWigMap := wig.MakeSkeleton(refChromSizes, 0)
	var curPosAlign string
	mapPostProbs := readPostProbs(postProbsFile)
	
	for idx := 0; idx < curLength; idx++ {
		curPosAlign = dna.BaseToString(mafInput[0].Species[0].SLine.Seq[idx]) + dna.BaseToString(mafInput[0].Species[1].SLine.Seq[idx]) + dna.BaseToString(mafInput[0].Species[2].SLine.Seq[idx]) + dna.BaseToString(mafInput[0].Species[3].SLine.Seq[idx])
		
		// hardcoded for 4 topologies
		outputWigMap["V0"].Values[idx] = mapPostProbs[curPosAlign][0]
		outputWigMap["V1"].Values[idx] = mapPostProbs[curPosAlign][1]
		outputWigMap["V2"].Values[idx] = mapPostProbs[curPosAlign][2]
		outputWigMap["V3"].Values[idx] = mapPostProbs[curPosAlign][3]
	}
	return outputWigMap
}