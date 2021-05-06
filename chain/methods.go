package chain

import (
	"github.com/vertgenlab/gonomics/fileio"
	"io"
)

//Perform operations on target, chains are one based

func (ch *Chain) GetChrom() string {
	return ch.TName
}

// Chains are 0-base, half open
func (ch *Chain) GetChromStart() int {
	if ch.TStrand {
		return ch.TStart
	} else {
		return getSwapTCoord(ch, true, false)
	}
}

func (ch *Chain) GetChromEnd() int {
	if ch.TStrand {
		return ch.TEnd
	} else {
		return getSwapTCoord(ch, false, true)
	}
}

/* Necessary function if chain ever needs to implement the Lift interface.
func (ch *Chain) UpdateLift(c string, start int, end int) {
	ch.TName = c
	if ch.TStrand {
		ch.TStart = start
	} else {

	}
	if ch.TStrand {
		ch.TEnd = end
	}
}*/

//WriteToFileHandle was added in order to implement the Interval and Lift interfaces.
func (ch *Chain) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, ch)
}

type ChainSlice []*Chain

func (ch ChainSlice) Len() int { return len(ch) }

func (ch ChainSlice) Swap(i, j int) { ch[i], ch[j] = ch[j], ch[i] }

func (ch *ChainSlice) Push(x interface{}) {
	answer := x.(*Chain)
	*ch = append(*ch, answer)
}

func (ch *ChainSlice) Pop() interface{} {
	oldQueue := *ch
	n := len(oldQueue)
	answer := oldQueue[n-1]
	*ch = oldQueue[:n-1]
	return answer
}

//TODO: Not sure how Dan wants this setup. What is the best way to implement the write in cases where we require headers and other arguments
/*
func (ch *Chain) WriteToFileHandle(file io.Writer, comments *HeaderComments) {
	WriteChain(file, ch, comments)
}

func (ch ChainSlice) Write(file string) {
	Write(file, ch)
}
*/

func (ch *Chain) NextRecord(file *fileio.EasyReader) bool {
	var done bool
	var next *Chain
	for next == nil && !done {
		next, done = NextChain(file)
	}
	if done {
		return true
	}
	*ch = *next
	return done
}

func (ch *Chain) Copy() interface{} {
	var answer *Chain = new(Chain)
	*answer = *ch
	return answer
}

type ByGenomicCoordinates struct {
	ChainSlice
}

func (g ByGenomicCoordinates) Less(i, j int) bool {
	// First sort criteria is chromosome
	if g.ChainSlice[i].GetChrom() < g.ChainSlice[j].GetChrom() {
		return true
	} else if g.ChainSlice[i].GetChrom() == g.ChainSlice[j].GetChrom() && g.ChainSlice[i].GetChromStart() < g.ChainSlice[j].GetChromStart() {
		// If chroms are equal then sort by start position
		return true
	} else if g.ChainSlice[i].GetChromStart() == g.ChainSlice[j].GetChromStart() && g.ChainSlice[i].GetChromEnd() < g.ChainSlice[j].GetChromEnd() {
		// If start positions are equal then the shorter region wins
		return true
	} else {
		return false
	}
}

func (ch *Chain) SwapBoth() *Chain {
	return SwapBoth(ch)
}
