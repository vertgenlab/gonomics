package chain

import (
	"io"
)

// GetChrom returns the name of the target chromosome.
func (ch Chain) GetChrom() string {
	return ch.TName
}

// GetChromStart returns the starting position on the target chromosome.
func (ch Chain) GetChromStart() int {
	if ch.TStrand {
		return ch.TStart
	} else {
		return getSwapTCoord(ch, true, false)
	}
}

// GetChromEnd returns the ending position on the target chromosome.
func (ch Chain) GetChromEnd() int {
	if ch.TStrand {
		return ch.TEnd
	} else {
		return getSwapTCoord(ch, false, true)
	}
}

// WriteToFileHandle writes the chain to the io.Writer
// This methods helps to implement the Interval and Lift interfaces.
func (ch Chain) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, ch)
}

// SwapBoth swaps the target and query and returns the chain after the swap.
func (ch Chain) SwapBoth() Chain {
	return SwapBoth(ch)
}
