package sam

import (
	"github.com/vertgenlab/gonomics/chromInfo"
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

	touched bool // true if Count or InsCount has been modified
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
	refmap := chromInfo.SliceToMap(header.Chroms)
	for read := range reads {
		if !passesFilters(read, filters) {
			continue
		}
		pb.sendPassed(read, header, includeNoData, refmap, send)
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

// sendPassed sends any Piles with a position before the start of s
func (pb *pileBuffer) sendPassed(s Sam, h Header, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile) {
	var done bool
	// starting from pb.idx to the end
	for i := pb.idx; i < len(pb.s); i++ {
		done = pb.checkSend(i, s, includeNoData, refmap, send)
		if done {
			return
		}
	}

	// starting from start to pb.idx
	for i := 0; i < pb.idx; i++ {
		done = pb.checkSend(i, s, includeNoData, refmap, send)
		if done {
			return
		}
	}
}

func (pb *pileBuffer) checkSend(i int, s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile) (done bool) {
	if pb.s[i].RefIdx == -1 {
		return
	}
	if pb.s[i].RefIdx < refmap[s.RName].Order {
		if pb.s[i].touched || includeNoData {
			send <- pb.s[i]
		}
		resetPile(&pb.s[i])
		return
	}
	if pb.s[i].Pos >= s.Pos {
		pb.idx = i
		return true
	}
	if pb.s[i].touched || includeNoData {
		send <- pb.s[i]
	}
	resetPile(&pb.s[i])
	return
}

// addBase to the pileBuffer at the input position
func (pb *pileBuffer) addBase(pos uint32, base dna.Base) {
	pb.getIdx(pos).Count[base]++
}

// addIns to the pileBuffer at the input position
func (pb *pileBuffer) addIns(pos uint32, seq string) {
	pb.getIdx(pos).InsCount[seq]++
}

// getIdx returns a pointer to the pile at index pos in s
func (pb *pileBuffer) getIdx(pos uint32) *Pile {
	queryIdx := int(pos - pb.s[pb.idx].Pos)
	if queryIdx > len(pb.s) {

	}
}

