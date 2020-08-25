package dnaThreeBit

import (
	"log"
	"math/bits"
)

func countRightMatches(seqOne []uint64, startOne int, lenOne int, seqTwo []uint64, startTwo int, lenTwo int) int {
	const bitsPerBase int = 3
	const bitsPerInt int = 64
	const basesPerInt int = bitsPerInt / bitsPerBase
	const ones uint64 = 0xFFFFFFFFFFFFFFFF

	var i, j int
	var seqDiff uint64 = 0
	var bitMatches, matches int = 0, 0

	var offsetOne int = (startOne % basesPerInt) * bitsPerBase
	var offsetTwo int = (startTwo % basesPerInt) * bitsPerBase

	if offsetOne != offsetTwo {
		log.Fatalf("Different offsets when comparing sequences\n")
	}

	i = startOne / basesPerInt
	j = startTwo / basesPerInt
	iEnd := (lenOne + basesPerInt - 1) / basesPerInt
	jEnd := (lenTwo + basesPerInt - 1) / basesPerInt

	seqDiff = seqOne[i] ^ seqTwo[j]
	seqDiff = seqDiff & (ones >> uint(offsetOne))

	bitMatches = bits.LeadingZeros64(seqDiff)
	matches = (bitMatches - offsetOne) / bitsPerBase

	for i, j = i+1, j+1; i < iEnd && j < jEnd && bitMatches == bitsPerInt; i, j = i+1, j+1 {
		seqDiff = seqOne[i] ^ seqTwo[j]
		bitMatches = bits.LeadingZeros64(seqDiff)
		matches += bitMatches / bitsPerBase
	}

	// TODO: I am not sure what to do here.  I am currently assuming that the "padding"
	// on the end of the two sequences is different so that it will not match
	//return common.Min(common.Min(matches, lenOne-startOne), lenTwo-startTwo)
	return matches
}

func countLeftMatches(seqOne []uint64, startOne int, seqTwo []uint64, startTwo int) int {

	const bitsPerBase int = 3
	const bitsPerInt int = 64
	const basesPerInt int = bitsPerInt / bitsPerBase
	const ones uint64 = 0xFFFFFFFFFFFFFFFF

	var i, j int
	var seqDiff uint64 = 0
	var bitMatches, matches int = 0, 0

	var offsetOne int = (startOne % basesPerInt) * bitsPerBase
	var offsetTwo int = (startTwo % basesPerInt) * bitsPerBase

	if offsetOne != offsetTwo {
		log.Fatalf("Different offsets when comparing sequences\n")
	}

	i = startOne / basesPerInt
	j = startTwo / basesPerInt

	seqDiff = seqOne[i] ^ seqTwo[j]
	seqDiff = seqDiff & (ones << uint(bitsPerInt-offsetOne-bitsPerBase))

	bitMatches = bits.TrailingZeros64(seqDiff)
	matches = (bitMatches - (bitsPerInt - offsetOne - bitsPerBase)) / bitsPerBase

	for i, j = i-1, j-1; i >= 0 && j >= 0 && bitMatches == bitsPerInt; i, j = i-1, j-1 {
		seqDiff = seqOne[i] ^ seqTwo[j]
		bitMatches = bits.TrailingZeros64(seqDiff)
		matches += bitMatches / bitsPerBase
	}

	return matches
}

