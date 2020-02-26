package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

type ThreeBit struct {
	Seq []uint64
	Len int
}

const (
	A          uint8 = 0
	C          uint8 = 1
	G          uint8 = 2
	T          uint8 = 3
	N          uint8 = 4
	PaddingOne uint8 = 5
	PaddingTwo uint8 = 6
)

func BasesToUint64(seq []dna.Base, start int, end int, padding uint8) uint64 {
	if end-start > 21 || start >= end {
		log.Fatalf("Error: when converting to ThreeBit. start=%d end=%d\n", start, end)
	}
	var idx int
	var answer uint64 = uint64(seq[start]) << 1
	for idx = start + 1; idx < end; idx++ {
		answer = answer << 3
		answer = answer | (uint64(seq[idx]) << 1)
	}
	for ; idx < start+21; idx++ {
		answer = answer << 3
		answer = answer | (uint64(padding) << 1)
	}
	return answer
}

func GetBase(frag *ThreeBit, pos uint) dna.Base {
	var lastBase uint64 = 7
	var idx uint = pos / 21
	var remainder uint = pos % 21
	var shift uint = 64 - 3*(remainder+1)
	return dna.Base((frag.Seq[idx] >> (shift)) & lastBase)
}

func NewThreeBit(inSeq []dna.Base, padding uint8) *ThreeBit {
	var sliceLenNeeded int = (len(inSeq) + 20) / 21
	var start, end int = 0, 0
	answer := ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)}
	for i := 0; i < sliceLenNeeded; i++ {
		start = i * 21
		end = common.Min(start+21, len(inSeq))
		answer.Seq[i] = BasesToUint64(inSeq, start, end, padding)
	}
	return (&answer)
}
