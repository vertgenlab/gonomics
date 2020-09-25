package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
)

func NewTwoBitRainbow(inSeq []dna.Base) []*TwoBit {
	var start, end, sliceLenNeeded int = 0, 0, 0
	clone := make([]dna.Base, len(inSeq))
	copy(clone, inSeq)
	answer := make([]*TwoBit, 32)
	for startOffset := 0; startOffset < 32; startOffset++ {
		sliceLenNeeded = (len(clone) + 31) / 32
		answer[startOffset] = &TwoBit{Seq: make([]uint64, sliceLenNeeded), Len: len(clone)}
		for i := 0; i < sliceLenNeeded; i++ {
			start = i * 32
			end = common.Min(start+32, len(clone))
			answer[startOffset].Seq[i] = BasesToUint64LeftAln(clone, start, end)
		}
		clone = append([]dna.Base{dna.A}, clone...)
	}
	return answer
}

// TwoBitRainbowDeReference is the same as NewTwoBitRainbow except it will be turn a non-pointer plice of dnaTwoBit
func TwoBitRainbowDeReference(inSeq []dna.Base) []TwoBit {
	var start, end, sliceLenNeeded int = 0, 0, 0
	clone := make([]dna.Base, len(inSeq))
	copy(clone, inSeq)
	answer := make([]TwoBit, 32)
	for startOffset := 0; startOffset < 32; startOffset++ {
		sliceLenNeeded = (len(clone) + 31) / 32
		answer[startOffset] = TwoBit{Seq: make([]uint64, sliceLenNeeded), Len: len(clone)}
		for i := 0; i < sliceLenNeeded; i++ {
			start = i * 32
			end = common.Min(start+32, len(clone))
			answer[startOffset].Seq[i] = BasesToUint64LeftAln(clone, start, end)
		}
		clone = append([]dna.Base{dna.A}, clone...)
	}
	return answer
}
