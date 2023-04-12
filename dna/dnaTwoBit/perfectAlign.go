package dnaTwoBit

import (
	"log"
	"math/bits"

	"github.com/vertgenlab/gonomics/numbers"
)

func CountRightMatches(one *TwoBit, startOne int, two *TwoBit, startTwo int) int {

	const bitsPerBase int = 2
	const bitsPerInt int = 64
	const basesPerInt int = bitsPerInt / bitsPerBase
	const ones uint64 = 0xFFFFFFFFFFFFFFFF

	var i, j int
	var seqDiff uint64 = 0
	var bitMatches, totalMatches int = 0, 0

	var offsetOne int = (startOne % basesPerInt) * bitsPerBase
	var offsetTwo int = (startTwo % basesPerInt) * bitsPerBase

	if offsetOne != offsetTwo {
		log.Fatalf("Error: Different offsets when comparing sequences\n")
	}

	i = startOne / basesPerInt
	j = startTwo / basesPerInt
	iEnd := (one.Len + basesPerInt - 1) / basesPerInt
	jEnd := (two.Len + basesPerInt - 1) / basesPerInt

	seqDiff = one.Seq[i] ^ two.Seq[j]
	seqDiff = seqDiff & (ones >> uint(offsetOne))

	bitMatches = bits.LeadingZeros64(seqDiff)
	totalMatches = bitMatches - offsetOne

	for i, j = i+1, j+1; i < iEnd && j < jEnd && bitMatches == bitsPerInt; i, j = i+1, j+1 {
		seqDiff = one.Seq[i] ^ two.Seq[j]
		bitMatches = bits.LeadingZeros64(seqDiff)
		totalMatches += bitMatches
	}

	// TODO: I am not sure what to do here.  The "padding" on the end of the sequence
	// when it does not fit nicely into 64 bits, is assumed to be different, so that
	// it will not match.  We could use the commented out Min() function instead
	return numbers.Min(numbers.Min(totalMatches/bitsPerBase, one.Len-startOne), two.Len-startTwo)
	//return matches
}

func CountLeftMatches(one *TwoBit, startOne int, two *TwoBit, startTwo int) int {
	const bitsPerBase int = 2
	const bitsPerInt int = 64
	const basesPerInt int = bitsPerInt / bitsPerBase
	const ones uint64 = 0xFFFFFFFFFFFFFFFF

	var i, j int
	var seqDiff uint64 = 0
	var bitMatches, totalMatches int = 0, 0

	var offsetOne int = (startOne % basesPerInt) * bitsPerBase
	var offsetTwo int = (startTwo % basesPerInt) * bitsPerBase

	if offsetOne != offsetTwo {
		log.Fatalf("Different offsets when comparing sequences\n")
	}
	var firstBitsNoLook int = bitsPerInt - offsetOne - bitsPerBase

	i = startOne / basesPerInt
	j = startTwo / basesPerInt

	seqDiff = one.Seq[i] ^ two.Seq[j]
	seqDiff = seqDiff & (ones << uint(firstBitsNoLook))

	bitMatches = bits.TrailingZeros64(seqDiff)
	totalMatches = bitMatches - (firstBitsNoLook)

	for i, j = i-1, j-1; i >= 0 && j >= 0 && bitMatches == bitsPerInt; i, j = i-1, j-1 {
		seqDiff = one.Seq[i] ^ two.Seq[j]
		bitMatches = bits.TrailingZeros64(seqDiff)
		totalMatches += bitMatches
	}

	return totalMatches / bitsPerBase
}
