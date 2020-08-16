package alleles

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"sync"
)
//TODO: Merge with countSam functions using interfaces???
func GoCountGirafAlleles(girafFilename string, reference *simpleGraph.SimpleGraph, minMapQ uint8) <-chan *GraphAllele {
	answer := make(chan *GraphAllele)
	var wg sync.WaitGroup
	wg.Add(1)
	go CountGirafAlleles(answer, girafFilename, reference, minMapQ, &wg)

	go func() {
		wg.Wait()
		close(answer)
	}()

	return answer
}

func CountGirafAlleles(answer chan<- *GraphAllele, girafFilename string, reference *simpleGraph.SimpleGraph, minMapQ uint8, wg *sync.WaitGroup) {
	girafChan := giraf.GoReadToChan(girafFilename)
	var currAlleles = make(map[GraphLocation]*AlleleCount)
	var runningCount = make([]*GraphLocation, 0)
	var progress int // TODO: Make option to print progress

	for read := range girafChan {
		runningCount = sendPassedPositionsGiraf(answer, read, girafFilename, runningCount, currAlleles)
		currAlleles, runningCount = countGraphRead(read, currAlleles, runningCount, reference, minMapQ, progress)
	}
	wg.Done()
}

func sendPassedPositionsGiraf(answer chan<- *GraphAllele, aln *giraf.Giraf, girafFilename string, runningCount []*GraphLocation, currAlleles map[GraphLocation]*AlleleCount) []*GraphLocation {
	for i := 0; i < len(runningCount); i++ {

		if runningCount[i].Node.Id < aln.Path.Nodes[0] {
			answer <- &GraphAllele{girafFilename, currAlleles[*runningCount[i]], runningCount[i]}
			delete(currAlleles, *runningCount[i])

			// Catch instance where every entry in running count is sent
			// Delete all of runningCount
			if i == len(runningCount)-1 {
				runningCount = nil
			}
			continue
		}

		if runningCount[i].Node.Id == aln.Path.Nodes[0] && int(runningCount[i].Pos) < (aln.Path.TStart) {
			answer <- &GraphAllele{girafFilename, currAlleles[*runningCount[i]], runningCount[i]}
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

func addPosToMap(currLocation *GraphLocation, currAlleles map[GraphLocation]*AlleleCount, runningCount []*GraphLocation) (map[GraphLocation]*AlleleCount, []*GraphLocation) {
	if _, ok := currAlleles[*currLocation]; !ok {
		currAlleles[*currLocation] = &AlleleCount{
			Ref: currLocation.Node.Seq[currLocation.Pos], Counts: 0, BaseAF: 0, BaseCF: 0, BaseGF: 0, BaseTF: 0, BaseAR: 0, BaseCR: 0, BaseGR: 0, BaseTR: 0, Indel: make([]Indel, 0)}
		runningCount = append(runningCount, currLocation)
	}
	return currAlleles, runningCount
}
