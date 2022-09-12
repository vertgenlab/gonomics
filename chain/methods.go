package chain

import (
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

//TODO: Not sure how Dan wants this setup. What is the best way to implement the write in cases where we require headers and other arguments
/*
func (ch *Chain) WriteToFileHandle(file io.Writer, comments *HeaderComments) {
	WriteChain(file, ch, comments)
}
*/

func (ch *Chain) SwapBoth() *Chain {
	return SwapBoth(ch)
}
