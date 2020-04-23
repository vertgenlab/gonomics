package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
)

// Syncs allele streams output from SamToAlleles and returns a []*Allele for each position
//TODO: still needs to tested with graph input
func SyncAlleleStreams(reference interface{}, alleleStreams ...<-chan *Allele) <-chan []*Allele {
	fmt.Println("\nSTART SYNC ALLELES")
	answer := make(chan []*Allele)

	go syncAlleles(reference, alleleStreams, answer)

	return answer
}

func syncAlleles(reference interface{}, alleleStreams []<-chan *Allele, answer chan<- []*Allele) {
	var i int
	var ok bool
	var refIdx int = 0
	var currLoc *Allele = &Allele{Location: &Location{"dummy", 0}}
	var currAlleles []*Allele = make([]*Allele, len(alleleStreams))
	var currAnswer []*Allele

	for nextPos(reference, currLoc, &refIdx) {
		for i = 0; i < len(alleleStreams); i++ {
			if currAlleles[i] == nil {
				currAlleles[i], ok = <-alleleStreams[i]
				if !ok {continue}
			}
			if *currAlleles[i].Location == *currLoc.Location {
				currAnswer = append(currAnswer, currAlleles[i])
				currAlleles[i] = nil
			}
		}

		if currAnswer != nil {
			answer <- currAnswer
		}
		currAnswer = nil
	}

	close(answer)
}

func nextPos (reference interface{}, currLoc *Allele, refIdx *int) bool {
	switch r := reference.(type) {
	case []*fasta.Fasta:

		if currLoc.Location.Chr == "dummy" {
			currLoc.Location.Chr = r[0].Name
			currLoc.Location.Pos = 0
			return true
		} else {
			if int(currLoc.Location.Pos) < len(r[*refIdx].Seq) { // Increment Position
				currLoc.Location.Pos++
				return true
			} else if *refIdx < (len(r) - 1) { // Increment Chromosome
				*refIdx++
				currLoc.Location.Chr = r[*refIdx].Name
				currLoc.Location.Pos = 0
				return true
			} else { // Report end of reference
				return false
			}
		}

	case *simpleGraph.SimpleGraph:

		if currLoc.Location.Chr == "dummy" {
			currLoc.Location.Chr = r.Nodes[0].Name
			currLoc.Location.Pos = 0
			return true
		} else {
			if int(currLoc.Location.Pos) < len(r.Nodes[*refIdx].Seq) { // Increment Position
				currLoc.Location.Pos++
				return true
			} else if *refIdx < len(r.Nodes) { // Increment Node
				*refIdx++
				currLoc.Location.Chr = r.Nodes[*refIdx].Name
				currLoc.Location.Pos = 0
				return true
			} else { // Report end of reference
				return false
			}
		}
	}
	return false
}