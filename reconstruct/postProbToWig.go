package reconstruct

import (
	"bufio"
	"encoding/csv"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/wig"
)

// readLines reads a whole file into memory and returns a slice of its lines
func ReadLines(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close() // TODO Riley said to hard close it, but it won't let me do that?
	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	file.Close()
	return lines, scanner.Err()
}

// readPostProbs reads the posterior probability CSV file output from TRAILS
// currently supports n_AB=n_ABC=3, which outputs 29 states
// and per state, sums over i for all V_#^i topology, for each n which gives four sums per state
// alternatively also supports a posterior probability file that only has 4 states
func readPostProbs(postProbsFilename string) map[string][]float64 {
	f, err := os.Open(postProbsFilename)
	if err != nil {
		log.Fatal("Unable to read input file "+postProbsFilename, err)
	}

	reader := csv.NewReader(f)
	mapPostProbs := make(map[string][]float64)
	topoSum := float64(0)

	uncondensedEndIdx := []int{10, 16, 22, 28}

	row, err := reader.Read()
	if err != nil {
		log.Fatal("Unable to read the first row from file", err)
	}

	condensedBool := false
	if len(row) == 5 {
		condensedBool = true
	} else if len(row) != 29 {
		log.Fatal("Can only take uncondensed or condensed posterior probability file for n_abc=3, n_ab=3.")
	}

	var val float64
	startIdx := 0
	if condensedBool {
		for {
			// hardcoded to ignore anything that isn't n_abc=3,n_ab=3 for 4 different topologies
			probs := make([]float64, 4)
			for idx := 1; idx <= 4; idx++ {
				probs[idx-1], err = strconv.ParseFloat(row[idx], 64)
				if err != nil {
					log.Fatalf("Couldn't convert float %s", row[idx])
				}
			}

			mapPostProbs[row[0]] = probs

			row, err = reader.Read()
			if err != nil {
				break
			}
		}
	} else {
		for {
			probs := make([]float64, 4)
			startIdx = 0
			for probIdx, endIdx := range uncondensedEndIdx {
				topoSum = 0
				for idx := startIdx + 2; idx <= endIdx; idx++ {
					val, err = strconv.ParseFloat(row[idx], 64)
					if err != nil {
						log.Fatalf("Couldn't convert float %s", row[idx])
					}
					topoSum += val
				}
				probs[probIdx] = topoSum
				startIdx = endIdx
			}

			mapPostProbs[row[1]] = probs

			row, err = reader.Read()
			if err != nil {
				break
			}
		}
	}

	f.Close()
	return mapPostProbs
}

// PostProbToWig converts a csv file of posterior probabilities for four topologies to a wig file with four tracks
// desiredSpecies MUST be in order, 1st listed species must be the MAF reference species
// outputs two wigs: outputWigMap is a wig containing the posterior probability for the whole sequence
// per topology; mapEstimateMap is a wig produced by comparing the same position across all topologies and
// taking the highest probability state as the MAP (maximum a posteriori) estimate
func PostProbToWig(postProbsFile string, mafInput []*maf.Maf, desiredSpecies []string) (map[string]wig.Wig, map[string]wig.Wig) {
	// make a chrom.sizes ChromInfo struct for the entire maf file
	if len(mafInput) == 0 {
		log.Fatal("mafInput is empty")
	}
	if len(desiredSpecies) == 0 {
		log.Fatal("desiredSpecies is empty")
	}
	if len(mafInput[0].Species) == 0 {
		log.Fatal("first maf block has no species")
	}
	if !strings.HasPrefix(mafInput[0].Species[0].SLine.Src, desiredSpecies[0]) {
		log.Fatal("Reference species ", desiredSpecies[0], " not found.")
	}

	refSize := mafInput[0].Species[0].SLine.SrcSize

	trackNames := []string{"V0", "V1", "V2", "V3"}
	chromSizes := make([]chromInfo.ChromInfo, 4)

	for idx, trackName := range trackNames {
		chromSizes[idx] = chromInfo.ChromInfo{
			Name:  trackName,
			Size:  refSize,
			Order: idx,
		}
	}

	outputWigMap := wig.MakeSkeleton(chromInfo.SliceToMap(chromSizes), 0)
	mapEstimateMap := wig.MakeSkeleton(chromInfo.SliceToMap(chromSizes), 0)

	// iterate through each maf block
	desiredSpeciesIdx := make([]int, 4)
	allSpeciesPresent := 0
	var curPosAlign string
	collectSpecies := make([]string, 0)
	mapPostProbs := readPostProbs(postProbsFile)
	validBases := regexp.MustCompile(`^[ACGTacgt]{4}$`)
	var adjustIdx int
	refSpeciesIdx := 0
	mapIdx := 0
	mapVal := float64(0)
	lookupAlignment := ""
	for _, block := range mafInput {
		allSpeciesPresent = 0
		collectSpecies = collectSpecies[:0]
		for desiredIdx, species := range desiredSpecies {
			for speciesIdx, _ := range block.Species {
				collectSpecies = append(collectSpecies, species)
				if strings.HasPrefix(block.Species[speciesIdx].Src, species) {
					desiredSpeciesIdx[desiredIdx] = speciesIdx
					allSpeciesPresent++
					break
				}
			}
		}
		if allSpeciesPresent == 4 {
			for idx := 0; idx < block.Species[refSpeciesIdx].SLine.Size; idx++ {
				curPosAlign = dna.BaseToString(block.Species[desiredSpeciesIdx[0]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[1]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[2]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[3]].SLine.Seq[idx])

				adjustIdx = block.Species[refSpeciesIdx].SLine.Start + idx
				if !validBases.MatchString(curPosAlign) {
					continue
				}

				// hardcoded for 4 topologies
				mapIdx = 0
				mapVal = 0
				lookupAlignment = strings.ToUpper(curPosAlign)
				alignmentProb, valid := mapPostProbs[lookupAlignment]
				if !valid || len(alignmentProb) != 4 {
					log.Fatal("No probability exists for ", lookupAlignment, ", check posterior probability file.")
				}

				for idx, val := range mapPostProbs[lookupAlignment] {
					outputWigMap[trackNames[idx]].Values[adjustIdx] = mapPostProbs[lookupAlignment][idx]
					// assign if state is the highest probability state and what is the probability
					if val > mapVal {
						mapVal = val
						mapIdx = idx
					}

				}
				mapEstimateMap[trackNames[mapIdx]].Values[adjustIdx] = mapPostProbs[lookupAlignment][mapIdx]
			}
		}
	}

	return outputWigMap, mapEstimateMap
}
