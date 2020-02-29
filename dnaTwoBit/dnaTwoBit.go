package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

type TwoBit struct {
	Seq []uint64
	Len int
}

const (
	A          uint8 = 0
	C          uint8 = 1
	G          uint8 = 2
	T          uint8 = 3
)

func BasesToUint64LeftAln(seq []dna.Base, start int, end int) uint64 {
	if end-start > 32 || end-start < 1 {
		log.Fatalf("Error: when converting to TwoBit. start=%d end=%d\n", start, end)
	}
	var answer uint64 = 0
	var i int = 0
	for ; i < end; i++ {
		answer = answer << 2
		answer = answer | uint64(seq[i])
	}
	for ; i < start+32; i++ {
		answer = answer << 2
	}
	return answer
}

func BasesToUint64RightAln(seq []dna.Base, start int, end int) uint64 {
        if end-start > 32 || end-start < 1 {
                log.Fatalf("Error: when converting to TwoBit. start=%d end=%d\n", start, end)
        }
        var answer uint64 = 0
        for ; start < end; start++ {
                answer = answer << 2
                answer = answer | uint64(seq[start])
        }
        return answer
}

func GetBase(frag *TwoBit, pos uint) dna.Base {
	var lastBase uint64 = 3
	var idx uint = pos / 32
	var remainder uint = pos % 32
	var shift uint = 64 - 2*(remainder+1)
	return dna.Base((frag.Seq[idx] >> (shift)) & lastBase)
}

func NewTwoBit(inSeq []dna.Base) *TwoBit {
	var sliceLenNeeded int = (len(inSeq) + 31) / 32
	var start, end int = 0, 0
	answer := TwoBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)}
	for i := 0; i < sliceLenNeeded; i++ {
		start = i * 32
		end = common.Min(start+32, len(inSeq))
		answer.Seq[i] = BasesToUint64LeftAln(inSeq, start, end)
	}
	return (&answer)
}

