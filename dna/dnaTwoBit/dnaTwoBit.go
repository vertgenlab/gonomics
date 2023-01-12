//Package dnaTwoBit implements data structures for two-bit encoding of DNA sequences.

package dnaTwoBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// TwoBit is a struct that encodes DNA sequences in two-bit compressed format.
// Seq contains the sequence information itself. Len specifies the length of the sequence.
type TwoBit struct {
	Seq []uint64
	Len int
}

const (
	A uint8 = 0
	C uint8 = 1
	G uint8 = 2
	T uint8 = 3
)

// BasesToUint64LeftAln converts a user-specified range of an input slice of dna.Base structs into a tw-bit
// encoded uint64 sequence. Sequences less than 32 bases long are aligned to the left of the uint64 memory space.
func BasesToUint64LeftAln(seq []dna.Base, start int, end int) uint64 {
	if end-start > 32 || end-start < 1 {
		log.Fatalf("Error: when converting to TwoBit. start=%d end=%d\n", start, end)
	}
	var answer uint64 = 0
	var i int = start
	for ; i < end; i++ {
		answer = answer << 2//left shift two positions
		answer = answer | uint64(seq[i])//bitwise OR, appends base to right end of answer.
	}
	for ; i < start+32; i++ {
		answer = answer << 2
	}
	return answer
}

// BasesToUint64RightAln converts a user-specified range of an input slice of dna.Base structs into a tw-bit
// encoded uint64 sequence. Sequences less than 32 bases long are aligned to the right of the uint64 memory space.
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

// GetBase decodes a dna.Base stuct from a position of an input TwoBit sequence.
func GetBase(frag *TwoBit, pos uint) dna.Base {
	var lastBase uint64 = 3
	var idx uint = pos / 32
	var remainder uint = pos % 32
	var shift uint = 64 - 2*(remainder+1)
	return dna.Base((frag.Seq[idx] >> (shift)) & lastBase)
}

// NewTwoBit parses a TwoBit sequence struct from an input slice of dna.Base structs. Remainder bases will be left-aligned.
func NewTwoBit(inSeq []dna.Base) *TwoBit {
	var sliceLenNeeded int = (len(inSeq) + 31) / 32
	var start, end int = 0, 0
	answer := TwoBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)}
	for i := 0; i < sliceLenNeeded; i++ {
		start = i * 32
		end = numbers.Min(start+32, len(inSeq))
		answer.Seq[i] = BasesToUint64LeftAln(inSeq, start, end)
	}
	return &answer
}
