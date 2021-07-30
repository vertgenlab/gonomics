package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// Pile stores the number of each base observed across multiple reads.
// The number of a particular base can be retrieved by using the dna.Base
// as the index for Pile.Count. E.g. the number of reads with A is Count[dna.A].
type Pile struct {
	RefIdx int
	Pos uint32
	Count [13]int // Count[dna.Base] == Number of observed dna.Base
	InsCount map[string]int // key is insertion sequence as string
}

// GoPileup inputs a channel of coordinate sorted Sam structs and generates
// a Pile for each base. The Pile is send through the return channel when
// the Pile position is no longer being updated.
//
// The input filter functions must all be true for a read to be included.
// All reads are included if filters == nil.
func GoPileup(reads <-chan Sam, header Header, includeNoData bool, filters []func(s Sam) bool) <-chan Pile {
	if header.Metadata.SortOrder[0] != Coordinate {
		log.Fatal("ERROR: (GoPileup) input sam/bam must be coordinate sorted")
	}
	pileChan := make(chan Pile, 1000)
	pileup(pileChan, reads, header, includeNoData, filters)
	return pileChan
}

// pileup generates multiple Pile based on the input reads and sends the Pile when passed.
func pileup(send chan<- Pile, reads <-chan Sam, header Header, includeNoData bool, filters []func(s Sam) bool) {
	pb := newPileBuffer()
	for read := range reads {
		if !passesFilters(read, filters) {
			continue
		}
		updatePile(pb, read)
	}
	close(send)
}

// passesFilters returns true if input Sam is true for all functions in filters.
func passesFilters(s Sam, filters []func(s Sam) bool) bool {
	for _, f := range filters {
		if !f(s) {
			return false
		}
	}
	return true
}

// updatePile updates the input pileBuffer with the data in s
func updatePile(pb *pileBuffer, s Sam) {
	var seqPos int
	refPos := s.Pos
	for i := range s.Cigar {
		switch s.Cigar[i].Op {
		case 'M', '=', 'X': // Match
			addMatch(pb, refPos, s.Seq[seqPos:seqPos + s.Cigar[i].RunLength])
			refPos++
			seqPos++

		case 'D': // Deletion
			addDeletion(pb, refPos, s.Cigar[i].RunLength)
			refPos += uint32(s.Cigar[i].RunLength)

		case 'I': // Insertion
			addInsertion(pb, refPos, s.Seq[seqPos:seqPos + s.Cigar[i].RunLength])
			seqPos += s.Cigar[i].RunLength

		default:
			if cigar.ConsumesReference(s.Cigar[i].Op) {
				refPos += uint32(s.Cigar[i].RunLength)
			}
			if cigar.ConsumesQuery(s.Cigar[i].Op) {
				seqPos += s.Cigar[i].RunLength
			}
		}
	}
}

// addMatch to pileBuffer
func addMatch(pb *pileBuffer, startPos uint32, seq []dna.Base) {
	var i uint32
	for i = range seq {
		pb.addBase(startPos + i, seq[i])
	}
}

// addDeletion to pileBuffer
func addDeletion(pb *pileBuffer, startPos uint32, length int) {
	var i uint32
	for i = 0; i < uint32(length); i++ {
		pb.addBase(startPos + i, dna.Gap)
	}
}

// addInsertion to pileBuffer
func addInsertion(pb *pileBuffer, startPos uint32, seq []dna.Base) {
	pb.addIns(startPos, dna.BasesToString(seq))
}


// ****** pileBuffer functions below ****** //

// pileBuffer stores a slice of piles for efficient reuse.
type pileBuffer struct {
	s []Pile
	idx int // index in s of lowest pos
}

func (pb *pileBuffer) addBase(pos uint32, base dna.Base) {

}

func (pb *pileBuffer) addIns(pos uint32, seq string) {

}

// newPileBuffer creates a newPileBuffer with len of 300
func newPileBuffer() *pileBuffer {
	pb := pileBuffer{
		s: make([]Pile, 300), // initialize to len of 2 standard reads
	}
	for i := range pb.s {
		pb.s[i].RefIdx = -1
	}
	return &pb
}

// resetPile sets a Pile to the default state
func resetPile(p *Pile) {
	p.RefIdx = -1
	p.Pos = 0
	for i := range p.Count {
		p.Count[i] = 0
	}
}

