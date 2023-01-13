package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// NewThreeBitRainbow builds a "rainbow table" of a sequence in ThreeBit format
// with every possible offset, so that there is always a version of the sequence
// that can be compared to another ThreeBit sequence using xor
func NewThreeBitRainbow(inSeq []dna.Base, padding ThreeBitBase) []*ThreeBit {
	var i, start, end, sliceLenNeeded, startOffset int
	answer := make([]*ThreeBit, 21)

	for startOffset = 0; startOffset < 21; startOffset++ {
		sliceLenNeeded = (len(inSeq) + startOffset + 20) / 21
		answer[startOffset] = &ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq) + startOffset}

		end = numbers.Min(21-startOffset, len(inSeq))
		answer[startOffset].Seq[0] = basesToUint64WithOffset(inSeq, 0, end, padding, startOffset)
		for i = 1; i < sliceLenNeeded; i++ {
			start = i*21 - startOffset
			end = numbers.Min(start+21, len(inSeq))
			answer[startOffset].Seq[i] = BasesToUint64(inSeq, start, end, padding)
		}
	}
	return answer
}
