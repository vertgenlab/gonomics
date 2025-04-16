package reconstruct

import (
	"github.com/vertgenlab/gonomics/wig"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"encoding/csv"
	// "fmt"
	"log"
	"strconv"
	"os"
	"strings"

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
		log.Fatal("PostProbtoWig expects only 1 maf block.")
	}

	// make a chrom.sizes ChromInfo struct for the entire maf file
	var chromSizeInfo []chromInfo.ChromInfo
	for idx, block := range mafInput {
		chromSizeInfo = append(chromSizeInfo, chromInfo.ChromInfo{Name: block.Species[0].Src, Size: block.Species[0].SLine.Size, order: idx})

	}
	chromSizes := SliceToMap(chromSizeInfo)
	outputWigMap := wig.MakeSkeleton(chromSizes, 0)

	// iterate through each maf block
	for _, block := range mafInput {
		//// iterate through the species in the 1st maf block and check that each species has the same sequence length
		// for idx, species := range block.Species {
		// 	curLength = species.SLine.Size
		// 	if idx == 0 {
		// 		prevLength = curLength
		// 	} else if prevLength != curLength {
		// 		log.Fatal("Species must have same sequence length.")
		// 	}
		// }
		
		//// TODO: ignore this block of code, want to take an input chrom sizes file
		//// make the chrom size file based on the maf input (all species should have the same chrom size)
		// chroms := make([]chromInfo.ChromInfo, 0)
		// var trackName string
		// for idx := 0; idx <= 3; idx++ {
		// 	trackName = fmt.Sprintf("%s%d", "V", idx)
		// 	chroms = append(chroms, chromInfo.ChromInfo{Name: trackName, Size: curLength, Order: idx})
		// }
		
		
		// read the posterior probability file and load the alignent at position -> probabilities for each topology mapping
		// refChromSizes := chromInfo.SliceToMap(chroms)
		var curPosAlign string
		mapPostProbs := readPostProbs(postProbsFile) // read the posterior probabilities for each of the 4 species 
													// (aligns each state with the 4 topologies)

		// iterate through the maf file and check the alignment for each of the four species
		// match it to the posterior state
		// and assign it to the correct idx in the wig, based on the reference species in the block 
		// (1st species, which should stay the same) for all the blocks
		// instead of doing the sequence length check, for each position, need to check that all four bases exist, 
		// otherwise leave as the default wig value (i.e. 0 across all tracks)
		for idx := block.Species.SLine.Start; idx < block.Species.SLine.Start + block.Species.SLine.Size; idx++ {
			curPosAlign = dna.BaseToString(block.Species[0].SLine.Seq[idx]) + dna.BaseToString(block.Species[1].SLine.Seq[idx]) + dna.BaseToString(block.Species[2].SLine.Seq[idx]) + dna.BaseToString(block.Species[3].SLine.Seq[idx])
			
			if strings.Contains(curPosAlign, "-") {
				continue
			}

			// hardcoded for 4 topologies
			outputWigMap["V0"].Values[idx] = mapPostProbs[curPosAlign][0]
			outputWigMap["V1"].Values[idx] = mapPostProbs[curPosAlign][1]
			outputWigMap["V2"].Values[idx] = mapPostProbs[curPosAlign][2]
			outputWigMap["V3"].Values[idx] = mapPostProbs[curPosAlign][3]
		}
	}
	
	return outputWigMap
}