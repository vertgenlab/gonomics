package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
)

// Pile stores the number of each base observed across multiple reads.
// The number of a particular base can be retrieved by using the dna.Base
// as the index for Pile.Count. E.g. the number of reads with A is Count[dna.A].
type Pile struct {
	RefIdx   int
	Pos      uint32
	Count    [13]int        // Count[dna.Base] == Number of observed dna.Base
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
	go pileup(pileChan, reads, header, includeNoData, filters)
	return pileChan
}

// pileup generates multiple Pile based on the input reads and sends the Pile when passed.
func pileup(send chan<- Pile, reads <-chan Sam, header Header, includeNoData bool, filters []func(s Sam) bool) {
	pb := newPileBuffer(300) // initialize to 2x std read size
	refmap := chromInfo.SliceToMap(header.Chroms)
	for read := range reads {
		if !passesFilters(read, filters) {
			continue
		}
		pb.sendPassed(read, includeNoData, refmap, send)
		updatePile(pb, read, refmap)
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
func updatePile(pb *pileBuffer, s Sam, refmap map[string]chromInfo.ChromInfo) {
	var seqPos int
	refPos := s.Pos
	refidx := refmap[s.RName].Order
	for i := range s.Cigar {
		switch s.Cigar[i].Op {
		case 'M', '=', 'X': // Match
			addMatch(pb, refidx, refPos, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength])
			refPos += uint32(s.Cigar[i].RunLength)
			seqPos += s.Cigar[i].RunLength

		case 'D': // Deletion
			addDeletion(pb, refidx, refPos, s.Cigar[i].RunLength)
			refPos += uint32(s.Cigar[i].RunLength)

		case 'I': // Insertion
			addInsertion(pb, refidx, refPos-1, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength])
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
func addMatch(pb *pileBuffer, refidx int, startPos uint32, seq []dna.Base) {
	for i := range seq {
		pb.addBase(refidx, startPos+uint32(i), seq[i])
	}
}

// addDeletion to pileBuffer
func addDeletion(pb *pileBuffer, refidx int, startPos uint32, length int) {
	var i uint32
	for i = 0; i < uint32(length); i++ {
		pb.addBase(refidx, startPos+i, dna.Gap)
	}
}

// addInsertion to pileBuffer
func addInsertion(pb *pileBuffer, refidx int, startPos uint32, seq []dna.Base) {
	pb.addIns(refidx, startPos, dna.BasesToString(seq))
}

// ****** pileBuffer functions below ****** //

// pileBuffer stores a slice of piles for efficient reuse.
type pileBuffer struct {
	s   []Pile // discontinuous positions are not allowed
	idx int    // index in s of lowest pos
}

// newPileBuffer creates a newPileBuffer with input size
func newPileBuffer(size int) *pileBuffer {
	pb := pileBuffer{
		s: make([]Pile, size), // initialize to len of 2 standard reads
	}
	for i := range pb.s {
		pb.s[i].RefIdx = -1
		pb.s[i].InsCount = make(map[string]int)
	}
	return &pb
}

// resetPile sets a Pile to the default state
func resetPile(p *Pile) {
	p.RefIdx = -1
	p.Pos = 0
	if !p.touched {
		return
	}
	for i := range p.Count {
		p.Count[i] = 0
	}
	p.InsCount = make(map[string]int)
	p.touched = false
}

// sendPassed sends any Piles with a position before the start of s
func (pb *pileBuffer) sendPassed(s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile) {
	if pb.s[pb.idx].RefIdx == -1 {
		return // buffer is empty
	}
	var done bool
	var lastPos uint32
	var lastRefIdx int
	// starting from pb.idx to the end
	for i := pb.idx; i < len(pb.s); i++ {
		lastPos = pb.s[i].Pos
		lastRefIdx = pb.s[i].RefIdx
		done = pb.checkSend(i, s, includeNoData, refmap, send)
		if done {
			return
		}
	}

	// starting from start to pb.idx
	for i := 0; i < pb.idx; i++ {
		lastPos = pb.s[i].Pos
		lastRefIdx = pb.s[i].RefIdx
		done = pb.checkSend(i, s, includeNoData, refmap, send)
		if done {
			return
		}
	}

	// if we have not returned by this point, there are positions with no data
	// between the last position send, and the beginning of s
	pb.idx = 0
	if !includeNoData {
		return
	}

	// send missing chromosome data
	dummyPile := Pile{}
	for lastRefIdx < refmap[s.RName].Order {
		for i := lastPos + 1; int(i) < refmap[s.RName].Size; i++ {
			dummyPile.RefIdx = lastRefIdx
			dummyPile.Pos = i
			send <- dummyPile
		}
		lastPos = 0
		lastRefIdx++
	}

	for lastRefIdx == refmap[s.RName].Order {
		for i := lastPos + 1; i < s.Pos; i++ {
			dummyPile.RefIdx = lastRefIdx
			dummyPile.Pos = i
			send <- dummyPile
		}
		lastRefIdx++
	}
}

// checkSend checks if pb.s[i] is before the start of s and sends if true. returns true when no more records to send.
func (pb *pileBuffer) checkSend(i int, s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile) (done bool) {
	if pb.s[i].RefIdx == -1 {
		return
	}
	if pb.s[i].RefIdx < refmap[s.RName].Order {
		if pb.s[i].touched || includeNoData {
			send <- pb.s[i]
		}
		pb.incrementIdx()
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
	pb.incrementIdx()
	resetPile(&pb.s[i])
	return
}

// addBase to the pileBuffer at the input position
func (pb *pileBuffer) addBase(refidx int, pos uint32, base dna.Base) {
	pile := pb.getIdx(refidx, pos)
	pile.Count[base]++
	pile.touched = true
}

// addIns to the pileBuffer at the input position
func (pb *pileBuffer) addIns(refidx int, pos uint32, seq string) {
	pile := pb.getIdx(refidx, pos)
	pile.InsCount[seq]++
	pile.touched = true
}

// getIdx returns a pointer to the pile at index pos in s
func (pb *pileBuffer) getIdx(refidx int, pos uint32) *Pile {
	if pos < pb.s[pb.idx].Pos {
		log.Panic("tried to retrieve past position. unsorted input?")
	}
	queryIdx := int(pos-pb.s[pb.idx].Pos) + pb.idx

	// calc if queryIdx wraps around to start of buffer
	if queryIdx >= len(pb.s) {
		queryIdx -= len(pb.s)
	}

	// check position for correct data, return if found
	if queryIdx < len(pb.s) && pb.s[queryIdx].Pos == pos {
		return &pb.s[queryIdx]
	}

	// at this point, one of the following is true:
	// 1: the buffer is empty and needs to be filled
	// 2: the buffer is partially full and can be filled to pos
	// 3: the buffer is full and needs to be expanded

	switch {
	case pb.s[pb.idx].RefIdx == -1: // buffer is empty (case 1)
		pb.initializeFromEmpty(pos, refidx)
		return &pb.s[pb.idx]

	case pb.s[pb.decrementIdx()].RefIdx == -1: // buffer is partially full (case 2)
		pb.fillFromPartial()
		return pb.getIdx(refidx, pos)

	default: // buffer is full (case 3)
		pb.expand()
		return pb.getIdx(refidx, pos)
	}

}

// initializeFromEmpty fills an empty pileBuffer starting with the position of s.Seq[0]
func (pb *pileBuffer) initializeFromEmpty(pos uint32, refidx int) {
	for ; pb.s[pb.idx].RefIdx == -1; pb.idx = pb.incrementIdx() {
		pb.s[pb.idx].Pos = pos
		pb.s[pb.idx].RefIdx = refidx
		pos++
	}
}

// fillFromPartial fills a partially filled pileBuffer
func (pb *pileBuffer) fillFromPartial() {
	start := pb.idx
	refidx := pb.s[pb.idx].RefIdx
	pos := pb.s[pb.idx].Pos + 1
	for pb.idx = pb.incrementIdx(); pb.idx != start; pb.idx = pb.incrementIdx() {
		if pb.s[pb.idx].RefIdx == -1 {
			pb.s[pb.idx].RefIdx = refidx
			pb.s[pb.idx].Pos = pos
		}
		pos++
	}
}

// expand the pileBuffer by 2x its current size
func (pb *pileBuffer) expand() {
	newpb := newPileBuffer(len(pb.s) * 2)
	for i := 0; i < len(pb.s); i++ {
		newpb.s[i] = pb.s[pb.idx]
		pb.idx = pb.incrementIdx()
	}
	*pb = *newpb
}

// incrementIdx returns the idx of the pile after pb.idx
func (pb *pileBuffer) incrementIdx() int {
	newidx := pb.idx + 1
	if newidx == len(pb.s) {
		return 0
	}
	return newidx
}

// decrementIdx returns the idx of the pile before pb.idx
func (pb *pileBuffer) decrementIdx() int {
	newidx := pb.idx - 1
	if newidx == -1 {
		return len(pb.s) - 1
	}
	return newidx
}

// String for debug
func (pb *pileBuffer) String() string {
	s := new(strings.Builder)
	s.WriteString(fmt.Sprintf(
		"\nCap: %d\n"+
			"Idx: %d\n"+
			"Data starting from 0:\n", len(pb.s), pb.idx))

	for i, val := range pb.s {
		s.WriteString(fmt.Sprintf("Idx: %d\t%s\n", i, val.String()))
	}

	return s.String()
}

// String for debug
func (p *Pile) String() string {
	return fmt.Sprintf("RefIdx: %d\tPos: %d\tCount: %v\tInsCount:%v", p.RefIdx, p.Pos, p.Count, p.InsCount)
}
