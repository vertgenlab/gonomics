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
	// touched is used as a quick check to see whether a pile struct
	// has been used. This value is used to avoid unnecessary allocations
	// in resetPile, and is used in checkSend to avoid sending zeroed structs
	// that happen to pass the user-defined filters.
}

// GoPileup inputs a channel of coordinate sorted Sam structs and generates
// a Pile for each base. The Pile is send through the return channel when
// the Pile position is no longer being updated.
//
// The input readFilter functions must all be true for a read to be included.
// The input pileFilter functions can be likewise be used to filter the output.
// All reads/piles are included if filters == nil.
func GoPileup(reads <-chan Sam, header Header, includeNoData bool, readFilters []func(s Sam) bool, pileFilters []func(p Pile) bool) <-chan Pile {
	if header.Metadata.SortOrder[0] != Coordinate {
		log.Fatal("ERROR: (GoPileup) input sam/bam must be coordinate sorted")
	}
	pileChan := make(chan Pile, 1000)
	go pileup(pileChan, reads, header, includeNoData, readFilters, pileFilters)
	return pileChan
}

// pileup generates multiple Pile based on the input reads and sends the Pile when passed.
func pileup(send chan<- Pile, reads <-chan Sam, header Header, includeNoData bool, readFilters []func(s Sam) bool, pileFilters []func(p Pile) bool) {
	pb := newPileBuffer(300) // initialize to 2x std read size
	refmap := chromInfo.SliceToMap(header.Chroms)
	for read := range reads {
		if !passesReadFilters(read, readFilters) {
			continue
		}
		pb.sendPassed(read, includeNoData, refmap, send, pileFilters)
		updatePile(pb, read, refmap)
	}
	close(send)
}

// passesReadFilters returns true if input Sam is true for all functions in filters.
func passesReadFilters(s Sam, filters []func(s Sam) bool) bool {
	for _, f := range filters {
		if !f(s) {
			return false
		}
	}
	return true
}

// passesPileFilters returns true if input Sam is true for all functions in filters.
func passesPileFilters(p Pile, filters []func(p Pile) bool) bool {
	for _, f := range filters {
		if !f(p) {
			return false
		}
	}
	return true
}

// updatePile updates the input pileBuffer with the data in buf
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
	buf []Pile // discontinuous positions are not allowed
	idx int    // index in buf of lowest pos
}

// newPileBuffer creates a newPileBuffer with input size
func newPileBuffer(size int) *pileBuffer {
	pb := pileBuffer{
		buf: make([]Pile, size), // initialize to len of 2 standard reads
	}
	for i := range pb.buf {
		pb.buf[i].RefIdx = -1
		pb.buf[i].InsCount = make(map[string]int)
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

// sendPassed sends any Piles with a position before the start of buf
func (pb *pileBuffer) sendPassed(s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile, pileFilters []func(p Pile) bool) {
	if pb.buf[pb.idx].RefIdx == -1 {
		return // buffer is empty
	}
	var done bool
	var lastPos uint32
	var lastRefIdx int
	// starting from pb.idx to the end
	for i := pb.idx; i < len(pb.buf); i++ {
		lastPos = pb.buf[i].Pos
		lastRefIdx = pb.buf[i].RefIdx
		done = pb.checkSend(i, s, includeNoData, refmap, send, pileFilters)
		if done {
			return
		}
	}

	// starting from start to pb.idx
	for i := 0; i < pb.idx; i++ {
		lastPos = pb.buf[i].Pos
		lastRefIdx = pb.buf[i].RefIdx
		done = pb.checkSend(i, s, includeNoData, refmap, send, pileFilters)
		if done {
			return
		}
	}

	// if we have not returned by this point, there are positions with no data
	// between the last position send, and the beginning of buf
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

// checkSend checks if pb.buf[i] is before the start of buf and sends if true. returns true when no more records to send.
func (pb *pileBuffer) checkSend(i int, s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile, pileFilters []func(p Pile) bool) (done bool) {
	if pb.buf[i].RefIdx == -1 {
		return
	}
	if pb.buf[i].RefIdx < refmap[s.RName].Order {
		if (pb.buf[i].touched && passesPileFilters(pb.buf[i], pileFilters)) || includeNoData {
			send <- pb.buf[i]
		}
		pb.incrementIdx()
		resetPile(&pb.buf[i])
		return
	}
	if pb.buf[i].Pos >= s.Pos {
		pb.idx = i
		return true
	}
	if (pb.buf[i].touched && passesPileFilters(pb.buf[i], pileFilters)) || includeNoData {
		send <- pb.buf[i]
	}
	pb.incrementIdx()
	resetPile(&pb.buf[i])
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

// getIdx returns a pointer to the pile at index pos in buf
func (pb *pileBuffer) getIdx(refidx int, pos uint32) *Pile {
	if pos < pb.buf[pb.idx].Pos {
		log.Panic("tried to retrieve past position. unsorted input?")
	}
	queryIdx := int(pos-pb.buf[pb.idx].Pos) + pb.idx

	// calc if queryIdx wraps around to start of buffer
	if queryIdx >= len(pb.buf) {
		queryIdx -= len(pb.buf)
	}

	// check position for correct data, return if found
	if queryIdx < len(pb.buf) && pb.buf[queryIdx].Pos == pos {
		return &pb.buf[queryIdx]
	}

	// at this point, one of the following is true:
	// 1: the buffer is empty and needs to be filled
	// 2: the buffer is partially full and can be filled to pos
	// 3: the buffer is full and needs to be expanded

	switch {
	case pb.buf[pb.idx].RefIdx == -1: // buffer is empty (case 1)
		pb.initializeFromEmpty(pos, refidx)
		return &pb.buf[pb.idx]

	case pb.buf[pb.decrementIdx()].RefIdx == -1: // buffer is partially full (case 2)
		pb.fillFromPartial()
		return pb.getIdx(refidx, pos)

	default: // buffer is full (case 3)
		pb.expand()
		return pb.getIdx(refidx, pos)
	}

}

// initializeFromEmpty fills an empty pileBuffer starting with the position of buf.Seq[0]
func (pb *pileBuffer) initializeFromEmpty(pos uint32, refidx int) {
	for ; pb.buf[pb.idx].RefIdx == -1; pb.idx = pb.incrementIdx() {
		pb.buf[pb.idx].Pos = pos
		pb.buf[pb.idx].RefIdx = refidx
		pos++
	}
}

// fillFromPartial fills a partially filled pileBuffer
func (pb *pileBuffer) fillFromPartial() {
	start := pb.idx
	refidx := pb.buf[pb.idx].RefIdx
	pos := pb.buf[pb.idx].Pos + 1
	for pb.idx = pb.incrementIdx(); pb.idx != start; pb.idx = pb.incrementIdx() {
		if pb.buf[pb.idx].RefIdx == -1 {
			pb.buf[pb.idx].RefIdx = refidx
			pb.buf[pb.idx].Pos = pos
		}
		pos++
	}
}

// expand the pileBuffer by 2x its current size
func (pb *pileBuffer) expand() {
	newpb := newPileBuffer(len(pb.buf) * 2)
	for i := 0; i < len(pb.buf); i++ {
		newpb.buf[i] = pb.buf[pb.idx]
		pb.idx = pb.incrementIdx()
	}
	*pb = *newpb
}

// incrementIdx returns the idx of the pile after pb.idx
func (pb *pileBuffer) incrementIdx() int {
	newidx := pb.idx + 1
	if newidx == len(pb.buf) {
		return 0
	}
	return newidx
}

// decrementIdx returns the idx of the pile before pb.idx
func (pb *pileBuffer) decrementIdx() int {
	newidx := pb.idx - 1
	if newidx == -1 {
		return len(pb.buf) - 1
	}
	return newidx
}

// String for debug
func (pb *pileBuffer) String() string {
	s := new(strings.Builder)
	s.WriteString(fmt.Sprintf(
		"\nCap: %d\n"+
			"Idx: %d\n"+
			"Data starting from 0:\n", len(pb.buf), pb.idx))

	for i, val := range pb.buf {
		s.WriteString(fmt.Sprintf("Idx: %d\t%buf\n", i, val.String()))
	}

	return s.String()
}

// String for debug
func (p *Pile) String() string {
	return fmt.Sprintf("RefIdx: %d\tPos: %d\tCount: %v\tInsCount:%v", p.RefIdx, p.Pos, p.Count, p.InsCount)
}
