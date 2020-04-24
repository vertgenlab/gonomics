package alleles

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/simpleGraph"
)

// Syncs allele streams output from SamToAlleles and returns a []*Allele with the Allele struct for each sample at position n
//TODO: still needs to tested with graph input
func SyncAlleleStreams(reference interface{}, memBufferSize int, alleleStreams ...<-chan *Allele) <-chan []*Allele {
	answer := make(chan []*Allele)
	go syncAlleles(reference, memBufferSize, alleleStreams, answer)
	return answer
}

func syncAlleles(reference interface{}, memBufferSize int, alleleStreams []<-chan *Allele, answer chan<- []*Allele) {
	var i int
	var refIdx int = 0
	// To get alleles for each position we increment throught the reference (whether linear or graph) and request alleles for each position
	// currLoc stores the current location we are asking for and begins with a dummy value to start for the nextPos function
	var currLoc *Allele = &Allele{Location: &Location{"dummy", 0}}
	// The currAlleles data structure is arranged as [SampleStream][Allele (length = memBufferSize)]
	// currAlleles transiently stores data read from each alleleStream and allows for a variable buffer size depending on the desired memory allocation
	var currAlleles [][]*Allele = make([][]*Allele, len(alleleStreams))
	var currAnswer []*Allele

	// Make the second dimension of currAlleles for each alleleStream
	for i = 0; i < len(alleleStreams); i++ {
		currAlleles[i] = make([]*Allele, 0)
	}

	// Loop through each position in the reference
	for nextPos(reference, currLoc, &refIdx) {
		// For each alleleStream
		for i = 0; i < len(alleleStreams); i++ {
			// If there is space in the memory buffer, read an allele from the stream and store it
			if len(currAlleles[i]) < memBufferSize {
				recdAllele, ok := <-alleleStreams[i]
				if !ok { // If reached EOF and the alleleStream has closed
					// AND all of the alleles from that stream have been sent
					if len(currAlleles[i]) == 0 {
						// Remove the closed alleleStream from the list and from currAlleles
						alleleStreams = append(alleleStreams[:i], alleleStreams[i+1:]...)
						currAlleles = append(currAlleles[:i], currAlleles[i+1:]...)
						// Decrement i since the current alleleStream corresponding to i has been deleted
						i--
					}
					continue
				} else {
					currAlleles[i] = append(currAlleles[i], recdAllele)
				}
			}

			// If the location we are asking for is the same as the first allele stored in currAlleles[i]
			if *currAlleles[i][0].Location == *currLoc.Location {
				// Append value to answer and delete allele from currAlleles
				currAnswer = append(currAnswer, currAlleles[i][0])
				currAlleles[i] = currAlleles[i][1:]
			}
		}

		if currAnswer != nil {
			answer <- currAnswer
		}
		currAnswer = nil
	}

	close(answer)
}

// Inputs a []*fasta.Fasta or a *simpleGraph.SimpleGraph and increments currLoc
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