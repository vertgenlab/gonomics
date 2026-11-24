package ilsReconstruct

import (
	"encoding/csv"
	"fmt"
	"log"
	"regexp"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/wig"
)

// readPostProbs reads the posterior probability CSV file output from TRAILS
// currently supports n_AB=n_ABC=3, which outputs 29 states
// and per state, sums over i for all V_#^i topology, for each n which gives four sums per state
// alternatively also supports a posterior probability file that only has 4 states
func readPostProbs(postProbsFilename string) map[string][]float64 {
	f := fileio.EasyOpen(postProbsFilename)

	reader := csv.NewReader(f)
	mapPostProbs := make(map[string][]float64)
	topoSum := float64(0)

	// Column groupings for an uncondensed posterior probability file with 29 TRAILS states
	// Creates groups of columns mapping to the four topologies (V0-3) for n_AB = n_ABC = 3.
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

// PostProbToWig converts TRAILS posterior probabilities into
// topology-specific WIG tracks using nucleotide patterns from a MAF
// multiple sequence alignment.
// desiredSpecies must contain N species in the same order used by
// TRAILS, with an ungapped reference species listed first. Will not
// give correct alignment if given gapped reference sequence.
// Returns two maps of WIG tracks: outputWigMap contains posterior
// probabilities for each topology, while mapEstimateMap contains
// the highest posterior probability assigned to its corresponding
// topology at each processed genomic position.
func PostProbToWig(postProbsFile string, mafInput []*maf.Maf, desiredSpecies []string, numTopos int) (map[string]wig.Wig, map[string]wig.Wig) {
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

	if numTopos <= 0 {
		log.Fatalf("numTopos must be a positive integer")
	}

	//// Use reference genome length to initialize four WIG tracks (one per topology)
	//// for posterior probabilities and MAP estimates.
	refSize := mafInput[0].Species[0].SLine.SrcSize

	chromSizes := make([]chromInfo.ChromInfo, numTopos)
	trackNames := make([]string, numTopos)
	for idx := 0; idx < numTopos; idx++ {
		trackNames[idx] = fmt.Sprintf("V%d", idx)

		chromSizes[idx] = chromInfo.ChromInfo{
			Name:  trackNames[idx],
			Size:  refSize,
			Order: idx,
		}
	}

	for idx, trackName := range trackNames {
		chromSizes[idx] = chromInfo.ChromInfo{
			Name:  trackName,
			Size:  refSize,
			Order: idx,
		}
	}

	// Initialize genome-length WIG tracks for each topolog
	// Unassigned positions retain posterior probability=zero
	// Want to storee all topology probabilities in outputWigMap, but record
	// only the highest-probability topology in mapEstimateMap
	outputWigMap := wig.MakeSkeleton(chromInfo.SliceToMap(chromSizes), 0)
	mapEstimateMap := wig.MakeSkeleton(chromInfo.SliceToMap(chromSizes), 0)

	// Locate desired species within each MAF block, as species
	// may appear in different orders or be absent from some blocks
	// Only process blocks containing all required species
	numSpecies := len(desiredSpecies)
	desiredSpeciesIdx := make([]int, numSpecies)

	allSpeciesPresent := 0
	var curPosAlign string
	collectSpecies := make([]string, 0)
	mapPostProbs := readPostProbs(postProbsFile)
	validBases := regexp.MustCompile(fmt.Sprintf(`^[ACGTacgt]{%d}$`, numSpecies))

	// adjustIdx determines reference genomic position corresponding to the
	// current alignment column for assigning posterior probabilities.
	var adjustIdx int
	refSpeciesIdx := 0
	mapIdx := 0
	mapVal := float64(0)
	lookupAlignment := ""
	for _, block := range mafInput {
		allSpeciesPresent = 0
		collectSpecies = collectSpecies[:0]
		for desiredIdx, species := range desiredSpecies {
			for speciesIdx := range block.Species {
				collectSpecies = append(collectSpecies, species)
				if strings.HasPrefix(block.Species[speciesIdx].Src, species) {
					desiredSpeciesIdx[desiredIdx] = speciesIdx
					allSpeciesPresent++
					break
				}
			}
		}

		if allSpeciesPresent == numSpecies {
			for idx := 0; idx < block.Species[refSpeciesIdx].SLine.Size; idx++ {
				adjustIdx = block.Species[refSpeciesIdx].SLine.Start + idx
				curPosAlign = dna.BaseToString(block.Species[desiredSpeciesIdx[0]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[1]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[2]].SLine.Seq[idx]) + dna.BaseToString(block.Species[desiredSpeciesIdx[3]].SLine.Seq[idx])

				if !validBases.MatchString(curPosAlign) {
					continue
				}

				mapIdx = 0
				mapVal = 0
				// Normalize nucleotide pattern to match keys in posterior probability lookup table
				lookupAlignment = strings.ToUpper(curPosAlign)
				alignmentProb, valid := mapPostProbs[lookupAlignment]
				if !valid || len(alignmentProb) != numTopos {
					log.Fatal("No probability exists for ", lookupAlignment, ", check posterior probability file.")
				}

				for idx, val := range mapPostProbs[lookupAlignment] {
					outputWigMap[trackNames[idx]].Values[adjustIdx] = mapPostProbs[lookupAlignment][idx]
					// assign if state is the highest probability state
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
