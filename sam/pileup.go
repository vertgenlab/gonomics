package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
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
	RefIdx    int
	Pos       uint32         // 1-base (like Sam)
	CountF    [13]int        // Count[dna.Base] == Number of observed dna.Base, forward reads
	CountR    [13]int        // Count[dna.Base] == Number of observed dna.Base
	InsCountF map[string]int // key is insertion sequence as string, value is number of observations, forward reads
	InsCountR map[string]int // key is insertion sequence as string, value is number of observations

	// Note that DelCountF/R DO NOT contribute to the total depth of a particular base.
	// They are only included to preserve multi-base deletion structure for downstream use.
	// Further note that DelCount is only recorded for the 5'-most base in the deletion.
	DelCountF map[int]int // key is the number of contiguous bases that are deleted, value is number of observations, forward reads
	DelCountR map[int]int // key is the number of contiguous bases that are deleted, value is number of observations, forward reads

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
func GoSyncPileups(samples ...<-chan Pile) <-chan []Pile {
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
			if buf[i].RefIdx != minRefIdx || buf[i].Pos != minPos || samples[i] == nil {
				data[i].RefIdx = -1 // mark as no data in output
				continue
			}

			data[i] = buf[i]              // save matching Pile to output data
			buf[i], isOpen = <-samples[i] // replace Pile in buf

			if !isOpen { // set channel to nil once closed
				buf[i].RefIdx = -1
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
		if (p[i].RefIdx >= 0 && p[i].RefIdx < p[minIdx].RefIdx) || (p[i].RefIdx == p[minIdx].RefIdx && p[i].Pos < p[minIdx].Pos) {
			minIdx = i
		}
	}
	//fmt.Println(minIdx, p)
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
	var lastSentRefIdx int
	var lastSentPos uint32
	start := newLinkedPileBuffer(300) // initialize to 2x std read size
	refmap := chromInfo.SliceToMap(header.Chroms)
	var read Sam
	for read = range reads {
		if !passesReadFilters(read, readFilters) {
			continue
		}
		sclipTerminalIns(&read)
		start, lastSentRefIdx, lastSentPos = sendPassedLinked(start, read, includeNoData, refmap, send, pileFilters, lastSentRefIdx, lastSentPos)
		updateLinkedPile(start, read, refmap)
	}
	sendRemaining(start, send, includeNoData, pileFilters, lastSentRefIdx, lastSentPos, refmap[read.RName])
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
	var readIsForward bool = !IsPaired(s) || IsForwardRead(s) // record unpaired as forward
	refPos := s.Pos
	refidx := refmap[s.RName].Order
	for i := range s.Cigar {
		switch s.Cigar[i].Op {
		case 'M', '=', 'X': // Match
			addMatchLinked(start, refidx, refPos, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength], readIsForward)
			refPos += uint32(s.Cigar[i].RunLength)
			seqPos += s.Cigar[i].RunLength

		case 'D': // Deletion
			addDeletionLinked(start, refidx, refPos, s.Cigar[i].RunLength, readIsForward)
			refPos += uint32(s.Cigar[i].RunLength)

		case 'I': // Insertion
			addInsertionLinked(start, refidx, refPos-1, s.Seq[seqPos:seqPos+s.Cigar[i].RunLength], readIsForward)
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

func addMatchLinked(start *Pile, refidx int, startPos uint32, seq []dna.Base, readIsForward bool) *Pile {
	for i := range seq {
		start = getPile(start, refidx, startPos+uint32(i))
		if readIsForward {
			start.CountF[seq[i]]++
		} else {
			start.CountR[seq[i]]++
		}
		start.touched = true
	}
	return start
}

func addDeletionLinked(start *Pile, refidx int, startPos uint32, length int, readIsForward bool) *Pile {
	start = getPile(start, refidx, startPos)
	if start.DelCountF == nil { // only make the map when we need it
		start.DelCountF = make(map[int]int)
	}
	if start.DelCountR == nil { // only make the map when we need it
		start.DelCountR = make(map[int]int)
	}
	if readIsForward {
		start.DelCountF[length]++
	} else {
		start.DelCountR[length]++
	}

	var i uint32
	for i = 0; i < uint32(length); i++ {
		start = getPile(start, refidx, startPos+i)
		if readIsForward {
			start.CountF[dna.Gap]++
		} else {
			start.CountR[dna.Gap]++
		}

		start.touched = true
	}
	return start
}

func addInsertionLinked(start *Pile, refidx int, startPos uint32, seq []dna.Base, readIsForward bool) *Pile {
	start = getPile(start, refidx, startPos)
	if start.InsCountF == nil { // only make the map when we need it
		start.InsCountF = make(map[string]int)
	}
	if start.InsCountR == nil { // only make the map when we need it
		start.InsCountR = make(map[string]int)
	}
	if readIsForward {
		start.InsCountF[dna.BasesToString(seq)]++
	} else {
		start.InsCountR[dna.BasesToString(seq)]++
	}
	start.touched = true
	return start
}

func getPile(start *Pile, refidx int, pos uint32) *Pile {
	if pos < start.Pos {
		log.Panicf("sent a record for writing before all the data was present. something went horribly wrong")
	}
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
			start = start.prev                 // back up to end of buffer
			expandLinkedPileBuffer(start, 300) // TODO: POTENTIAL MEMORY LEAK. HARD CAP???
			start = start.next                 // move into start of newly added buffer

		default:
			start = start.next
		}
	}
	return start
}

func sendRemaining(start *Pile, send chan<- Pile, includeNoData bool, pileFilters []func(p Pile) bool, lastRefIdx int, lastPos uint32, ref chromInfo.ChromInfo) {
	for start.RefIdx != -1 {
		if (start.touched || includeNoData) && passesPileFilters(*start, pileFilters) {
			send <- *start
		}
		start.RefIdx = -1
		start = start.next
	}

	// if we have not returned by this point, there are positions with no data
	// between the last position send, and the beginning of buf
	if !includeNoData {
		return
	}

	// send missing chromosome data
	dummyPile := Pile{}
	for i := lastPos + 1; i <= uint32(ref.Size); i++ {
		dummyPile.RefIdx = lastRefIdx
		dummyPile.Pos = i
		send <- dummyPile
	}
}

func sendPassedLinked(start *Pile, s Sam, includeNoData bool, refmap map[string]chromInfo.ChromInfo, send chan<- Pile, pileFilters []func(p Pile) bool, lastRefIdx int, lastPos uint32) (newStart *Pile, lastSentRefIdx int, lastSentPos uint32) {
	for start.RefIdx != refmap[s.RName].Order || start.Pos < s.Pos-1 { // the -1 on s.Pos is to handle cases where a read begins with an insertion
		if start.RefIdx == -1 {
			break
		}
		if (start.touched || includeNoData) && passesPileFilters(*start, pileFilters) {
			lastRefIdx = start.RefIdx
			lastPos = start.Pos

			send <- *start
		}

		resetPile(start)
		start = start.next
	}

	// if we have not returned by this point, there are positions with no data
	// between the last position send, and the beginning of buf
	if !includeNoData {
		return start, lastRefIdx, lastPos
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

	return start, lastRefIdx, lastPos
}

// resetPile sets a Pile to the default state
func resetPile(p *Pile) {
	p.RefIdx = -1
	p.Pos = 0
	if !p.touched {
		return
	}
	for i := range p.CountF {
		p.CountF[i] = 0
		p.CountR[i] = 0
	}

	p.InsCountF = nil // only make maps when we need them
	p.InsCountR = nil
	p.DelCountF = nil
	p.DelCountR = nil
	p.touched = false
}

// sclipTerminalIns will convert an insertion on the left or right end of the read to a soft clip
func sclipTerminalIns(s *Sam) {
	if len(s.Cigar) == 0 || s.Cigar[0].Op == '*' {
		return
	}
	if s.Cigar[0].Op == 'I' {
		s.Cigar[0].Op = 'S'
	}
	if s.Cigar[len(s.Cigar)-1].Op == 'I' {
		s.Cigar[len(s.Cigar)-1].Op = 'S'
	}
}

// String for debug
func (p *Pile) String() string {
	return fmt.Sprintf("RefIdx: %d\tPos: %d\tCountF: %v\tCountR: %v\tInsCountF: %v\tInsCountR: %v\tDelCountF: %v\tDelCountR: %v\tNext: %p\tPrev: %p",
		p.RefIdx, p.Pos, p.CountF, p.CountR, p.InsCountF, p.InsCountR, p.DelCountF, p.DelCountR, p.next, p.prev)
}

// consensusType represents the possible types of consensus variants in a Pile,
// which includes Bases, INDELs, or undefined (in the case where the Pile contains no data).
type consensusType byte

const (
	Base      consensusType = 0
	Insertion consensusType = 1
	Deletion  consensusType = 2
	Undefined consensusType = 3
)

// Consensus contains information on the most observed variant in a Pile struct.
type Consensus struct {
	Base      dna.Base
	Insertion []dna.Base
	Deletion  int
	Type      consensusType
}

// PileConsensus returns a Consensus struct from an input Pile struct
// representing information about the consensus variant observed in that Pile.
func PileConsensus(p Pile) Consensus {
	// first we check if consensus is a base
	var max int = p.CountF[dna.A] + p.CountR[dna.A]
	tiedConsensus := make([]Consensus, 1)
	tiedConsensus[0] = Consensus{Base: dna.A, Type: Base}

	max, tiedConsensus = getMaxBase(p, max, dna.C, tiedConsensus)
	max, tiedConsensus = getMaxBase(p, max, dna.G, tiedConsensus)
	max, tiedConsensus = getMaxBase(p, max, dna.T, tiedConsensus)
	//question for Dan: main version checks dna.Gap here, is that necessary if we check the deletion field?
	max, tiedConsensus = getMaxInsertion(p, max, tiedConsensus)
	max, tiedConsensus = getMaxDeletion(p, max, tiedConsensus)

	if max < 1 {
		return Consensus{Type: Undefined}
	}
	if len(tiedConsensus) == 1 {
		return tiedConsensus[0]
	}
	return tiedConsensus[numbers.RandIntInRange(0, len(tiedConsensus))]
}

// helper function of PileConsensus, finds the max deletion in a Pile struct.
func getMaxDeletion(p Pile, currMax int, tiedConsensus []Consensus) (int, []Consensus) {
	var count, i int
	var inMap bool
	var seenOnPosStrand = make(map[int]int, 0)
	for i = range p.DelCountF {
		seenOnPosStrand[i] = 1
		count = p.DelCountF[i]
		if _, inMap = p.DelCountR[i]; inMap {
			count += p.DelCountR[i]
		}
		if count > currMax {
			tiedConsensus = tiedConsensus[:1]
			tiedConsensus[0] = Consensus{Deletion: i, Type: Deletion}
			currMax = count
		}
		if count == currMax {
			tiedConsensus = append(tiedConsensus, Consensus{Deletion: i, Type: Deletion})
		}
	}

	for i = range p.DelCountR {
		if _, inMap = seenOnPosStrand[i]; !inMap {
			count = p.DelCountR[i]
			if count > currMax {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0] = Consensus{Deletion: i, Type: Deletion}
				currMax = count
			}
			if count == currMax {
				tiedConsensus = append(tiedConsensus, Consensus{Deletion: i, Type: Deletion})
			}
		}
	}
	return currMax, tiedConsensus
}

// helper function of PileConsensus, finds the max insertion in a Pile struct.
func getMaxInsertion(p Pile, currMax int, tiedConsensus []Consensus) (int, []Consensus) {
	var count int
	var inMap bool
	var seenOnPosStrand = make(map[string]int, 0)
	for i := range p.InsCountF {
		seenOnPosStrand[i] = 1
		count = p.InsCountF[i]
		if _, inMap = p.InsCountR[i]; inMap {
			count += p.InsCountR[i]
		}
		if count > currMax {
			tiedConsensus = tiedConsensus[:1]
			tiedConsensus[0] = Consensus{Insertion: dna.StringToBases(i), Type: Insertion}
			currMax = count
		}
		if count == currMax {
			tiedConsensus = append(tiedConsensus, Consensus{Insertion: dna.StringToBases(i), Type: Insertion})
		}
	}

	//we have to check the p.InsCountR map for any insertions observed only on the minus strand.
	for i := range p.InsCountR {
		if _, inMap = seenOnPosStrand[i]; !inMap {
			count = p.InsCountR[i] //we don't have to sum since we are guaranteed not to have seen this insertion on the positive strand.
			if count > currMax {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0] = Consensus{Insertion: dna.StringToBases(i), Type: Insertion}
				currMax = count
			}
			if count == currMax {
				tiedConsensus = append(tiedConsensus, Consensus{Insertion: dna.StringToBases(i), Type: Insertion})
			}
		}
	}

	return currMax, tiedConsensus
}

// getMaxBase is a helper function of PileConsensus, finds the consensus base in a Pile struct.
func getMaxBase(p Pile, currMax int, testBase dna.Base, tiedConsensus []Consensus) (int, []Consensus) {
	var count int = p.CountF[testBase] + p.CountR[testBase]
	if count > currMax {
		tiedConsensus = tiedConsensus[:1] // reset tied bases
		tiedConsensus[0] = Consensus{Base: testBase, Type: Base}
		return count, tiedConsensus
	}

	if count == currMax {
		tiedConsensus = append(tiedConsensus, Consensus{Base: testBase, Type: Base})
	}

	return currMax, tiedConsensus
}
