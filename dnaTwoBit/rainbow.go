package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// TwoBitRainbow takes in a dna sequences and returns the "rainbow"
// of that sequence in TwoBit format.  This consists of all possible
// shifts of the original sequence so that when this sequence must
// becompared to another seq, a bit-wise xor can be used.
func TwoBitRainbow(inSeq []dna.Base) []TwoBit {
	var start, end, sliceLenNeeded int = 0, 0, 0
	clone := make([]dna.Base, len(inSeq))
	copy(clone, inSeq)
	answer := make([]TwoBit, 32)
	for startOffset := 0; startOffset < 32; startOffset++ {
		sliceLenNeeded = (len(clone) + 31) / 32
		answer[startOffset] = TwoBit{Seq: make([]uint64, sliceLenNeeded), Len: len(clone)}
		for i := 0; i < sliceLenNeeded; i++ {
			start = i * 32
			end = numbers.Min(start+32, len(clone))
			answer[startOffset].Seq[i] = BasesToUint64LeftAln(clone, start, end)
		}
		clone = append([]dna.Base{dna.A}, clone...)
	}
	return answer
}
