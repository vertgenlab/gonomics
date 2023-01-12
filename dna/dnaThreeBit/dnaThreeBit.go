package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

// ThreeBit is a struct to represent long DNA sequences in a memory efficient format
// Seq holds the encoded DNA sequence.  The left-most three bits of Seq[0] hold the
// first base.  The first base is *not* encoded in the three least significant bits.
// Because 64 is not divisible by three, the right-most (least significant) bit
// is not used in the encoding scheme.  Len is the length of the DNA sequence, not
// the length of the Seq slice.
type ThreeBit struct {
	Seq []uint64
	Len int
}

// This could have been uint8, but with that I found lots
// of casting between uint64 and uint8.  I don't think many
// ThreeBitBases will be sitting around for a long time by themselves
// so I don't think the extra memory will be noticed.
// Even though it is encoded as 64 bits, only the last three can be used (zero to seven)
type ThreeBitBase uint64

// The four bases {A,C,G,T} and N are encoded with the numbers 0 to 4.
// If the length of the DNA sequence stored in the threeBit is not divisible by
// 21, then there will be "left over" space in Seq[len(Seq)-1].  Currently,
// if that "left over" space is filled with PaddingOne for one threeBit,
// and PaddingTwo for another threeBit, the threeBits can be quickly compared
// for perfect matches.  The padding is so that the bits in those bases will not match.
// TODO: the padding implementation is a bit messy and should probably be changed in
// the future.  This could be done by having the perfect match functions compare the
// number of matches they will return to the theoretical maximum based on the sequence
// lengths or padding our a mask could be applied in the comparison function.
const (
	A          ThreeBitBase = 0
	C          ThreeBitBase = 1
	G          ThreeBitBase = 2
	T          ThreeBitBase = 3
	N          ThreeBitBase = 4
	PaddingOne ThreeBitBase = 5
	PaddingTwo ThreeBitBase = 6
)

// BasesToUint64 will take a section of seq from start to end (left-closed, right open)
// that is not more than 21 bases
// and return the uint64 encoding of it that would be used in a ThreeBit.
// padding can be either PaddingOne or PaddingTwo.  If sequences will be compared
// to each other, they should have different padding values.
func BasesToUint64(seq []dna.Base, start int, end int, padding ThreeBitBase) uint64 {
	if end-start > 21 || start >= end {
		log.Fatalf("Error: when converting to ThreeBit. start=%d end=%d\n", start, end)
	}
	var idx int
	var answer uint64 = 0
	for idx = start; idx < end; idx++ {
		answer = answer << 3 // not needed first time through, but does not hurt
		answer = answer | (uint64(seq[idx]) << 1)
	}
	for ; idx < start+21; idx++ {
		answer = answer << 3
		answer = answer | (uint64(padding) << 1)
	}
	return answer
}

// basesToUint64WithOffset places padding at the beginning of the sequence (as well as the end, which is normal)
// so that sequences may be "in register" with each other so that an xor can quickly compare them for equality
func basesToUint64WithOffset(seq []dna.Base, start int, end int, padding ThreeBitBase, offset int) uint64 {
	if end-start+offset > 21 || start >= end {
		log.Fatalf("Error: when converting to ThreeBit. start=%d end=%d\n", start, end)
	}
	var idx int
	var answer uint64 = 0
	for idx = 0; idx < offset; idx++ {
		answer = answer << 3 // not needed the first time through the loop, but does not hurt
		answer = answer | (uint64(padding) << 1)
	}
	for idx = start; idx < end; idx++ {
		answer = answer << 3
		answer = answer | (uint64(seq[idx]) << 1)
	}
	for ; idx < start+21-offset; idx++ {
		answer = answer << 3
		answer = answer | (uint64(padding) << 1)
	}
	return answer
}

// GetThreeBitBase returns a value equivalent to a single base within the ThreeBit at position, pos.
func GetThreeBitBase(fragment *ThreeBit, pos int) ThreeBitBase {
	if pos < 0 || pos >= fragment.Len {
		log.Fatalf("Error: asked for base at position:%d for a sequence with length:%d\n", pos, fragment.Len)
	}
	var upos = uint(pos)
	const lastBase uint64 = 7 // right-most three bits are one
	var idx uint = upos / 21
	var remainder uint = upos % 21
	var shift uint = 64 - 3*(remainder+1)
	return ThreeBitBase((fragment.Seq[idx] >> shift) & lastBase)
}

// GetBase returns a value equivalent to a single base within the ThreeBit at position, pos.
func GetBase(fragment *ThreeBit, pos int) dna.Base {
	return dna.Base(GetThreeBitBase(fragment, pos))
}

// NewThreeBit creates a ThreeBit encoding of inSeq with padding on the end
func NewThreeBit(inSeq []dna.Base, padding ThreeBitBase) *ThreeBit {
	var sliceLenNeeded int = (len(inSeq) + 20) / 21
	var start, end int = 0, 0
	answer := ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq)}
	for i := 0; i < sliceLenNeeded; i++ {
		start = i * 21
		end = numbers.Min(start+21, len(inSeq))
		answer.Seq[i] = BasesToUint64(inSeq, start, end, padding)
	}
	return &answer
}

// newThreeBitWithOffset creates a ThreeBit encoding of inSeq with an offset of "offset" positions
// at the start of the sequence so that sequences may be in register with each other.
func newThreeBitWithOffset(inSeq []dna.Base, padding ThreeBitBase, offset int) *ThreeBit {
	var sliceLenNeeded int = (len(inSeq) + offset + 20) / 21
	var start, end int = 0, 0
	answer := ThreeBit{Seq: make([]uint64, sliceLenNeeded), Len: len(inSeq) + offset}
	answer.Seq[0] = basesToUint64WithOffset(inSeq, 0, numbers.Min(len(inSeq), 21-offset), padding, offset)
	for i := 1; i < sliceLenNeeded; i++ {
		start = i*21 - offset
		end = numbers.Min(start+21, len(inSeq))
		answer.Seq[i] = BasesToUint64(inSeq, start, end, padding)
	}
	return &answer
}
