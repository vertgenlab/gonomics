package alleles

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"sync"
)

// sendPassedPositions sends positions that have been passed in the file
// TODO: Version for graph. Maybe readin graph and use simpleGraph.GetSortOrder to determine ordering
func sendPassedPositionsSam(answer chan<- *Allele, aln *sam.SamAln, samFilename string, runningCount []*Location, currAlleles map[Location]*AlleleCount) []*Location {
	for i := 0; i < len(runningCount); i++ {

		if runningCount[i].Chr != aln.RName {
			answer <- &Allele{samFilename, currAlleles[*runningCount[i]], runningCount[i]}
			delete(currAlleles, *runningCount[i])

			// Catch instance where every entry in running count is sent
			// Delete all of runningCount
			if i == len(runningCount)-1 {
				runningCount = nil
			}
			continue
		}

		if runningCount[i].Pos < (aln.Pos - 1) {
			answer <- &Allele{samFilename, currAlleles[*runningCount[i]], runningCount[i]}
			delete(currAlleles, *runningCount[i])

			// Catch instance where every entry in running count is sent
			// Delete all of runningCount
			if i == len(runningCount)-1 {
				runningCount = nil
			}

		} else {
			// Remove sent values from count
			runningCount = runningCount[i:]
			break
		}
	}
	return runningCount
}

func countSamRead(aln *sam.SamAln, currAlleles map[Location]*AlleleCount, runningCount []*Location, ref map[string][]dna.Base, minMapQ int64, progress int) (map[Location]*AlleleCount, []*Location) {

	if aln.Cigar[0].Op == '*' {
		return currAlleles, runningCount
	}

	if aln.MapQ < minMapQ {
		return currAlleles, runningCount
	}

	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var i, j, k int
	var currentIndel Indel
	var indelSeq []dna.Base
	var OrigRefIndex int64
	var Match bool

	// Count the bases
	progress++
	SeqIndex = 0
	RefIndex = aln.Pos - 1

	for i = 0; i < len(aln.Cigar); i++ {
		currentSeq = aln.Seq

		if aln.Cigar[i].Op == 'D' {
			OrigRefIndex = RefIndex
			indelSeq = make([]dna.Base, 1)

			// First base in indel is the base prior to the indel sequence per VCF standard format
			indelSeq[0] = ref[aln.RName][OrigRefIndex-1]

			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {

				// If the position has already been added to the map, move along
				_, ok := currAlleles[Location{aln.RName, RefIndex}]

				// If the position is NOT in the map, initialize
				if !ok {
					currAlleles[Location{aln.RName, RefIndex}] = &AlleleCount{
						Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
					runningCount = append(runningCount, &Location{aln.RName, RefIndex})
				}

				// Keep track of deleted sequence
				indelSeq = append(indelSeq, ref[aln.RName][RefIndex])

				currAlleles[Location{aln.RName, RefIndex}].Counts++
				RefIndex++
			}

			Match = false
			for j = 0; j < len(currAlleles[Location{aln.RName, OrigRefIndex}].Indel); j++ {
				// If the deletion has already been seen before, increment the existing entry
				// For a deletion the indelSeq should match the Ref
				if dna.CompareSeqsIgnoreCase(indelSeq, currAlleles[Location{aln.RName, OrigRefIndex}].Indel[j].Ref) == 0 &&
					dna.CompareSeqsIgnoreCase(indelSeq[:1], currAlleles[Location{aln.RName, OrigRefIndex}].Indel[j].Alt) == 0 {
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, OrigRefIndex}].Indel[j].CountF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, OrigRefIndex}].Indel[j].CountR++
					}

					Match = true
					break
				}
			}

			// If the deletion has not been seen before, then append it to the Del slice
			// For Alt indelSeq[:1] is used to give me a slice of just the first base in the slice which we defined earlier
			if Match == false {

				currentIndel = Indel{indelSeq, indelSeq[:1], 0, 0}
				if sam.IsForwardRead(aln) == true {
					currentIndel.CountF++
				} else if sam.IsReverseRead(aln) == false {
					currentIndel.CountR++
				}
				currAlleles[Location{aln.RName, OrigRefIndex}].Indel = append(currAlleles[Location{aln.RName, OrigRefIndex}].Indel, currentIndel)
			}

			//Handle insertion relative to ref
			//The base after the inserted sequence is annotated with an Ins read
		} else if aln.Cigar[i].Op == 'I' {

			// If the position has already been added to the map, move along
			_, ok := currAlleles[Location{aln.RName, RefIndex}]

			// If the position is NOT in the map, initialize
			if !ok {
				currAlleles[Location{aln.RName, RefIndex}] = &AlleleCount{
					Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
				runningCount = append(runningCount, &Location{aln.RName, RefIndex})
			}

			// Loop through read sequence and keep track of the inserted bases
			indelSeq = make([]dna.Base, 1)

			// First base in indel is the base prior to the indel sequence per VCF standard format
			indelSeq[0] = ref[aln.RName][RefIndex-1]

			// Keep track of inserted sequence by moving along the read
			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {
				indelSeq = append(indelSeq, currentSeq[SeqIndex])
				SeqIndex++
			}

			Match = false
			for j = 0; j < len(currAlleles[Location{aln.RName, RefIndex}].Indel); j++ {
				// If the inserted sequence matches a previously inserted sequence, then increment the count
				// For an insertion, the indelSeq should match the Alt
				if dna.CompareSeqsIgnoreCase(indelSeq, currAlleles[Location{aln.RName, RefIndex}].Indel[j].Alt) == 0 &&
					dna.CompareSeqsIgnoreCase(indelSeq[:1], currAlleles[Location{aln.RName, RefIndex}].Indel[j].Ref) == 0 {
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].Indel[j].CountF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].Indel[j].CountR++
					}
					Match = true
					break
				}
			}

			if Match == false {
				currentIndel = Indel{indelSeq[:1], indelSeq, 0, 0}
				if sam.IsForwardRead(aln) == true {
					currentIndel.CountF++
				} else if sam.IsReverseRead(aln) == true {
					currentIndel.CountR++
				}
				currAlleles[Location{aln.RName, RefIndex}].Indel = append(currAlleles[Location{aln.RName, RefIndex}].Indel, currentIndel)
			}

			// Note: Insertions do not contribute to the total counts as the insertion is associated with the previous reference base

			//Handle matching pos relative to ref
		} else if cigar.CigarConsumesReference(*aln.Cigar[i]) {

			for k = 0; k < int(aln.Cigar[i].RunLength); k++ {

				//if the position has already been added to the matrix, move along
				_, ok := currAlleles[Location{aln.RName, RefIndex}]

				//if the position is NOT in the matrix, add it
				if !ok {
					currAlleles[Location{aln.RName, RefIndex}] = &AlleleCount{
						Ref: ref[aln.RName][RefIndex], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
					runningCount = append(runningCount, &Location{aln.RName, RefIndex})
				}

				switch currentSeq[SeqIndex] {
				case dna.A:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseAF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseAR++
					}
					currAlleles[Location{aln.RName, RefIndex}].Counts++
				case dna.T:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseTF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseTR++
					}
					currAlleles[Location{aln.RName, RefIndex}].Counts++
				case dna.G:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseGF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseGR++
					}
					currAlleles[Location{aln.RName, RefIndex}].Counts++
				case dna.C:
					if sam.IsForwardRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseCF++
					} else if sam.IsReverseRead(aln) == true {
						currAlleles[Location{aln.RName, RefIndex}].BaseCR++
					}
					currAlleles[Location{aln.RName, RefIndex}].Counts++
				}
				SeqIndex++
				RefIndex++
			}
		} else if aln.Cigar[i].Op != 'H' {
			SeqIndex = SeqIndex + aln.Cigar[i].RunLength
		}
	}
	return currAlleles, runningCount
}

func GoCountSamAlleles(samFilename string, reference []*fasta.Fasta, minMapQ int64) <-chan *Allele {
	answer := make(chan *Allele)
	var wg sync.WaitGroup
	wg.Add(1)
	go CountSamAlleles(answer, samFilename, reference, minMapQ, &wg)

	go func() {
		wg.Wait()
		close(answer)
	}()

	return answer
}

func CountSamAlleles(answer chan<- *Allele, samFilename string, reference []*fasta.Fasta, minMapQ int64, wg *sync.WaitGroup) {
	samChan, _ := sam.GoReadToChan(samFilename)
	var currAlleles = make(map[Location]*AlleleCount)
	var runningCount = make([]*Location, 0)
	var progress int // TODO: Make option to print progress

	fasta.AllToUpper(reference)
	ref := fasta.FastaMap(reference)

	for read := range samChan {
		runningCount = sendPassedPositionsSam(answer, read, samFilename, runningCount, currAlleles)
		currAlleles, runningCount = countSamRead(read, currAlleles, runningCount, ref, minMapQ, progress)
	}
	wg.Done()
}

func GoCountGirafAlleles(girafFilename string, reference *simpleGraph.SimpleGraph, minMapQ uint8) <-chan *Allele {
	answer := make(chan *Allele)
	var wg sync.WaitGroup
	wg.Add(1)

	//TODO: simpleGraph.AllToUpper(reference)

	go CountGirafAlleles(answer, girafFilename, reference, minMapQ, &wg)

	go func() {
		wg.Wait()
		close(answer)
	}()

	return answer
}

func CountGirafAlleles(answer chan<- *Allele, girafFilename string, reference *simpleGraph.SimpleGraph, minMapQ uint8, wg *sync.WaitGroup) {
	girafChan := giraf.GoReadToChan(girafFilename)
	var currAlleles = make(map[GraphLocation]*AlleleCount)
	var runningCount = make([]*GraphLocation, 0)
	var progress int // TODO: Make option to print progress

	for read := range girafChan {
		runningCount = sendPassedPositionsGiraf(runningCount)
		currAlleles, runningCount = countGraphRead(read, currAlleles, runningCount, reference, minMapQ, progress)
	}
	wg.Done()
}

func sendPassedPositionsGiraf(runningCount []*GraphLocation) []*GraphLocation {

	return runningCount
}

func countBase(base dna.Base, indel *Indel, aln *giraf.Giraf, currLoc *GraphLocation, currAlleles map[GraphLocation]*AlleleCount) map[GraphLocation]*AlleleCount {
	currAlleles[*currLoc].Counts++

	if indel == nil { // If counting an SNV
		switch base {
		case dna.A:
			if giraf.IsReverseRead(aln) {
				currAlleles[*currLoc].BaseAR++
			} else { // If reads are unpaired then all counts will be recorded as fwd
				currAlleles[*currLoc].BaseAF++
			}
		case dna.C:
			if giraf.IsReverseRead(aln) {
				currAlleles[*currLoc].BaseCR++
			} else {
				currAlleles[*currLoc].BaseCF++
			}
		case dna.G:
			if giraf.IsReverseRead(aln) {
				currAlleles[*currLoc].BaseGR++
			} else {
				currAlleles[*currLoc].BaseGF++
			}
		case dna.T:
			if giraf.IsReverseRead(aln) {
				currAlleles[*currLoc].BaseTR++
			} else {
				currAlleles[*currLoc].BaseTF++
			}
		}
	} else { // If counting an indel
		for j := 0; j < len(currAlleles[*currLoc].Indel); j++ {
			// If the indel has already been seen before, increment the existing entry
			if dna.CompareSeqsIgnoreCase(indel.Ref, currAlleles[*currLoc].Indel[j].Ref) == 0 &&
				dna.CompareSeqsIgnoreCase(indel.Alt, currAlleles[*currLoc].Indel[j].Alt) == 0 {
				if giraf.IsForwardRead(aln) == true {
					currAlleles[*currLoc].Indel[j].CountF++
				} else if giraf.IsReverseRead(aln) == true {
					currAlleles[*currLoc].Indel[j].CountR++
				}
				return currAlleles
			}
		}

		// If the indel has not been seen before, then append it to the Indel slice
		if giraf.IsForwardRead(aln) == true {
			indel.CountF++
		} else if giraf.IsReverseRead(aln) == false {
			indel.CountR++
		}
		currAlleles[*currLoc].Indel = append(currAlleles[*currLoc].Indel, *indel)
	}

	return currAlleles
}

func countGraphRead(aln *giraf.Giraf, currAlleles map[GraphLocation]*AlleleCount, runningCount []*GraphLocation, ref *simpleGraph.SimpleGraph, minMapQ uint8, progress int) (map[GraphLocation]*AlleleCount, []*GraphLocation) {
	if aln.Aln[0].Op == '*' {
		return currAlleles, runningCount
	}

	if aln.MapQ < minMapQ {
		return currAlleles, runningCount
	}

	var i, refPos int64
	var readNodeIdx int = 0
	var currNode *simpleGraph.Node = ref.Nodes[aln.Path.Nodes[readNodeIdx]]
	var currLoc, originalLoc GraphLocation
	var currIndel Indel

	for _, cig := range aln.Aln {
		switch cig.Op {
		case '=': // Match
			for i = 0; i < cig.RunLength; i++ { // For each consecutive matching base
				currLoc = GraphLocation{Node: currNode, Pos: refPos}                               // define location on reference
				currAlleles, runningCount = addPosToMap(&currLoc, currAlleles, runningCount)       // add the pos to map if necessary
				currAlleles = countBase(currAlleles[currLoc].Ref, nil, aln, &currLoc, currAlleles) // count the base
				refPos++                                                                           // move to the next base in the ref
				if int(refPos) > len(currNode.Seq)-1 {                                             // check if next base is on another node, if so move to next node
					readNodeIdx++
					currNode = ref.Nodes[aln.Path.Nodes[readNodeIdx]]
					refPos = 0
				}
			}

		case 'X': // Mismatch
			for i = 0; i < cig.RunLength; i++ { // For each consecutive mismatched base
				currLoc = GraphLocation{Node: currNode, Pos: refPos}                         // define location on reference
				currAlleles, runningCount = addPosToMap(&currLoc, currAlleles, runningCount) // add the pos to map if necessary
				currAlleles = countBase(cig.Sequence[i], nil, aln, &currLoc, currAlleles)    // count the base
				refPos++                                                                     // move to the next base in the ref
				if int(refPos) > len(currNode.Seq)-1 {                                       // check if next base is on another node, if so move to next node
					readNodeIdx++
					currNode = ref.Nodes[aln.Path.Nodes[readNodeIdx]]
					refPos = 0
				}
			}

		case 'I': // Insertion
			// First base in indel is the base prior to the indel sequence per VCF standard format
			currLoc = GraphLocation{Node: currNode, Pos: refPos - 1}                     // define location on reference
			currAlleles, runningCount = addPosToMap(&currLoc, currAlleles, runningCount) // add the pos to map if necessary

			currIndel = Indel{ // initialize indel entry
				Ref:    make([]dna.Base, 1),
				Alt:    make([]dna.Base, 1),
				CountF: 0,
				CountR: 0,
			}

			currIndel.Ref[0] = currLoc.Node.Seq[currLoc.Pos] // add the base prior the the deleted sequence per vcf standard
			currIndel.Alt[0] = currLoc.Node.Seq[currLoc.Pos]

			currIndel.Alt = append(currIndel.Alt, cig.Sequence...) // add inserted bases to the Alt field

			currAlleles = countBase(0, &currIndel, aln, &currLoc, currAlleles) // count indel to the base prior to insertion. Base argument will be ignored in this call

		case 'D': // Deletion
			// First base in indel is the base prior to the indel sequence per VCF standard format
			currLoc = GraphLocation{Node: currNode, Pos: refPos - 1}                     // define location on reference
			originalLoc = currLoc                                                        // save the start location for use later (in case we jump nodes adding deleted sequence)
			currAlleles, runningCount = addPosToMap(&currLoc, currAlleles, runningCount) // add the pos to map if necessary

			currIndel = Indel{ // initialize indel entry
				Ref:    make([]dna.Base, 1),
				Alt:    make([]dna.Base, 1),
				CountF: 0,
				CountR: 0,
			}
			currIndel.Ref[0] = currLoc.Node.Seq[currLoc.Pos] // add the base prior the the deleted sequence per vcf standard
			currIndel.Alt[0] = currLoc.Node.Seq[currLoc.Pos]

			for i = 0; i < cig.RunLength; i++ { // add deleted sequence to currIndel.Ref
				//TODO: Think carefully about whether a read count should be added for each deleted base or if only one count should be recorded at the base prior to deletion
				if int(refPos) > len(currNode.Seq)-1 { // check if next base is on another node, if so move to next node
					readNodeIdx++
					currNode = ref.Nodes[aln.Path.Nodes[readNodeIdx]]
					refPos = 0
				}
				currIndel.Ref = append(currIndel.Ref, currLoc.Node.Seq[refPos]) // add base to currIndel
				refPos++                                                        // move to next base
			}
			currAlleles = countBase(0, &currIndel, aln, &originalLoc, currAlleles) // count indel to the base prior to deletion (originalLoc). Base argument will be ignored in this call

		case 'M': // Likely that cigar is using sam format and not fancy format. Throw error to let user know they may need to reformat.
			log.Fatalf("ERROR: Unrecognized character 'M' in cigar. Reformat to fancy cigars to use this functionality.")
		default: //TODO handle other cigar runes, maybe with cigar.ConsumesReference/Query
			log.Fatalf("ERROR: Unrecognized character '%v' in cigar for read %s. Valid characters are '=' 'X' 'I' 'D'.", cig.Op, aln.QName)
		}
	}
	return currAlleles, runningCount
}

func addPosToMap(currLocation *GraphLocation, currAlleles map[GraphLocation]*AlleleCount, runningCount []*GraphLocation) (map[GraphLocation]*AlleleCount, []*GraphLocation) {
	if _, ok := currAlleles[*currLocation]; !ok {
		currAlleles[*currLocation] = &AlleleCount{
			Ref: currLocation.Node.Seq[currLocation.Pos], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
		runningCount = append(runningCount, currLocation)
	}
	return currAlleles, runningCount
}
