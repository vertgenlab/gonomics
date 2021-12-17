package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// Pile stores the number of each base observed across multiple reads.
// The number of a particular base can be retrieved by using the dna.Base
// as the index for Pile.Count. E.g. the number of reads with A is Count[dna.A].
//
// Linked list mechanics:
// In the pileup function the Pile struct is part of a circular doubly-linked
// list in which Pile.next/prev is the next/prev contiguous position. Discontiguous
// positions are not allowed. When adding base information to the pile struct, the
// pointer to the start of the read position should be kept, and a second pointer
// will increment along the list as bases are added. The pointer to the beginning
// of the list should only be incremented when a position is passed and set to zero.
type Pile struct {
	RefIdx   int
	Pos      uint32
	Count    [13]int        // Count[dna.Base] == Number of observed dna.Base
	InsCount map[string]int // key is insertion sequence as string, value is number of observations

	touched bool // true if Count or InsCount has been modified
	// touched is used as a quick check to see whether a pile struct
	// has been used. This value is used to avoid unnecessary allocations
	// in resetPile, and is used in checkSend to avoid sending zeroed structs
	// that happen to pass the user-defined filters.

	next *Pile
	prev *Pile
}

// GoSyncPileups inputs a slice of channels receiving Pile structs and syncs
// their outputs to a new output channel that organizes Piles into a slice.
// The returned slice maintains the order of the input slice. Any of the input
// samples that have no data for the returned position will have the RefIdx
// field set to -1.
func GoSyncPileups(samples []<-chan Pile) <-chan []Pile {
	synced := make(chan []Pile, 1000)
	go syncPileups(samples, synced)
	return synced
}

// syncPileups performs syncs positions from an input slice of pile channels
// to an output channel of []Pile.
func syncPileups(samples []<-chan Pile, output chan<- []Pile) {
	var allClosed, isOpen bool
	buf := make([]Pile, len(samples))

	// fill buffer with 1 element from each channel
	for i := range samples {
		buf[i], isOpen = <-samples[i]
		if !isOpen {
			buf[i].RefIdx = -1
			samples[i] = nil
		}
	}

	// loop until all channels are closed
	var minRefIdx int
	var minPos uint32
	for !allClosed {
		data := make([]Pile, len(samples))
		minRefIdx, minPos = getMinPileCoords(buf)
		for i := range buf { // check for matching records in the buffer
			if buf[i].RefIdx != minRefIdx || buf[i].Pos != minPos {
				data[i].RefIdx = -1 // mark as no data in output
				continue
			}
			data[i] = buf[i]              // save matching Pile to output data
			buf[i], isOpen = <-samples[i] // replace Pile in buf

			if !isOpen { // set channel to nil once closed
				samples[i] = nil
				if allChannelsClosed(samples) { // each time a channel is closed, check if they are all closed
					allClosed = true
				}
			}
		}
		output <- data
	}
	close(output)
}

// getMinPileCoords returns the RefIdx and Pos of the Pile with the lowest coordinates.
func getMinPileCoords(p []Pile) (int, uint32) {
	var minIdx int
	for i := 1; i < len(p); i++ {
		if p[minIdx].RefIdx == -1 { // handles case where first index has no data
			minIdx = i
			continue
		}
		if p[i].RefIdx >= 0 && p[i].RefIdx <= p[minIdx].RefIdx && p[i].Pos < p[minIdx].Pos {
			minIdx = i
		}
	}
	return p[minIdx].RefIdx, p[minIdx].Pos
}

// allChanelsClosed returns true if all input channels are set to nil (closed).
func allChannelsClosed(c []<-chan Pile) bool {
	for i := range c {
		if c[i] != nil {
			return false
		}
	}
	return true
}

// GoPileup inputs a channel of coordinate sorted Sam structs and generates
// a Pile for each base. The Pile is sent through the return channel when
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
	go pileupLinked(pileChan, reads, header, includeNoData, readFilters, pileFilters)
	return pileChan
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

func pileupLinked(send chan<- Pile, reads <-chan Sam, header Header, includeNoData bool, readFilters []func(s Sam) bool, pileFilters []func(p Pile) bool) {
	start := newLinkedPileBuffer(300) // initialize to 2x std read size
	refmap := chromInfo.SliceToMap(header.Chroms)
	for read := range reads {
		if !passesReadFilters(read, readFilters) {
			continue
		}
		start = sendPassedLinked(start, read, includeNoData, refmap, send, pileFilters)
		updateLinkedPile(start, read, refmap)
	}
	close(send)
}

// newLinkedPilebuffer creates size Pile structs as a doubly-linked list and returns the start.
func newLinkedPileBuffer(size int) (start *Pile) {
	buf := make([]Pile, size)

	for i := range buf {
		buf[i].RefIdx = -1

		switch i {
		case 0: // first in slice
			buf[i].next = &buf[i+1]
			buf[i].prev = &buf[len(buf)-1]

		case len(buf) - 1: // last in slice
			buf[i].next = &buf[0]
			buf[i].prev = &buf[i-1]

		default:
			buf[i].next = &buf[i+1]
			buf[i].prev = &buf[i-1]
		}
	}
	return &buf[0]
}

// expandLinkedPileBuffer adds toAdd new Piles to the list directly after start.
func expandLinkedPileBuffer(start *Pile, toAdd int) {
	end := start.next
	newStart := newLinkedPileBuffer(toAdd)
	newEnd := newStart.prev

	start.next = newStart
	newStart.prev = start
	end.prev = newEnd
	newEnd.next = end
}

// updateLinkedPile updates the linked list with the data in s.
func updateLinkedPile(start *Pile, s Sam, refmap map[string]chromInfo.ChromInfo) {
	var seqPos int
	refPos := s.Pos
	refidx := refmap[s.RName].Order
	for i := range s.Cigar {
		switch s.Cigar[i].Op {
		case 'M', '=', 'X': // Match
			addMatchLinked(start, refidx, refPos, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength])
			refPos += uint32(s.Cigar[i].RunLength)
			seqPos += s.Cigar[i].RunLength

		case 'D': // Deletion
			addDeletionLinked(start, refidx, refPos, s.Cigar[i].RunLength)
			refPos += uint32(s.Cigar[i].RunLength)

		case 'I': // Insertion
			addInsertionLinked(start, refidx, refPos-1, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength])
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

func addMatchLinked(start *Pile, refidx int, startPos uint32, seq []dna.Base) *Pile {
	for i := range seq {
		start = getPile(start, refidx, startPos+uint32(i))
		start.Count[seq[i]]++
		start.touched = true
	}
	return start
}

func addDeletionLinked(start *Pile, refidx int, startPos uint32, length int) *Pile {
	var i uint32
	for i = 0; i < uint32(length); i++ {
		start = getPile(start, refidx, startPos+i)
		start.Count[dna.Gap]++
		start.touched = true
	}
	return start
}

func addInsertionLinked(start *Pile, refidx int, startPos uint32, seq []dna.Base) *Pile {
	start = getPile(start, refidx, startPos)
	if start.InsCount == nil { // only make the map when we need it
		start.InsCount = make(map[string]int)
	}
	start.InsCount[dna.BasesToString(seq)]++
	start.touched = true
	return start
}

func getPile(start *Pile, refidx int, pos uint32) *Pile {
	for start.RefIdx != refidx || start.Pos != pos {
		switch {
		case start.prev.RefIdx == -1 && start.RefIdx == -1: // no data in buffer, start at pos
			start.RefIdx = refidx
			start.Pos = pos
			return start

		case start.RefIdx == -1: // previous has data, current is empty. update from previous
			start.RefIdx = start.prev.RefIdx
			start.Pos = start.prev.Pos + 1

		case start.Pos < start.prev.Pos: // looped to start of buffer. need to expand
			start = start.prev // back up to end of buffer
			expandLinkedPileBuffer(start, 300)
			start = start.next // move into start of newly added buffer

		default:
			start = start.next
		}
	}
	return start
}

func sendPassedLinked(start *Pile, s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile, pileFilters []func(p Pile) bool) (newStart *Pile) {
	var lastRefIdx int
	var lastPos uint32
	for start.RefIdx != refmap[s.RName].Order || start.Pos != s.Pos {
		if start.RefIdx == -1 {
			break
		}

		if (start.touched || includeNoData) && passesPileFilters(*start, pileFilters) {
			lastRefIdx = start.RefIdx
			lastPos = start.Pos

			send <- *start
			resetPile(start)
			start = start.next
		}
	}

	// if we have not returned by this point, there are positions with no data
	// between the last position send, and the beginning of buf
	if !includeNoData {
		return start
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

	return start
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

	p.InsCount = nil // only make map when we need it
	p.touched = false
}

// String for debug
func (p *Pile) String() string {
	return fmt.Sprintf("RefIdx: %d\tPos: %d\tCount: %v\tInsCount:%v\tNext:%p\tPrev:%p", p.RefIdx, p.Pos, p.Count, p.InsCount, p.next, p.prev)
}
