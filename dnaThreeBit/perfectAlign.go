package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/common"
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

	// TODO: I am not sure what to do here.  The "padding" on the end of the sequence
	// when it does not fit nicely into 64 bits, is assumed to be different, so that
	// it will not match.  We could use the commented out Min() function instead
	return common.Min(common.Min(matches, lenOne-startOne), lenTwo-startTwo)
	//return matches
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

/*
func countRightMatches(seqOne []uint64, startOne uint, lenOne uint, seqTwo []uint64, startTwo uint, lenTwo uint) int {

	const bitsPerBase uint = 3
	const bitsPerInt = 64
	const basesPerInt = bitsPerInt / bitsPerBase

	var offsetOne uint = (startOne % basesPerInt) * bitsPerBase
	var offsetTwo uint = (startTwo % basesPerInt) * bitsPerBase

	if offsetTwo < offsetOne {
		seqOne, seqTwo = seqTwo, seqOne
		startOne, startTwo = startTwo, startOne
		lenOne, lenTwo = lenTwo, lenOne
		offsetOne, offsetTwo = offsetTwo, offsetOne
	}

	var rightShift uint = offsetTwo - offsetOne
	var leftShift uint = bitsPerInt - rightShift
	var i, j uint
	var seqDiff, prevBits uint64 = 0, 0
	const ones uint64 = 0xFFFFFFFFFFFFFFFF
	var matches int = 0

	iStart := startOne/basesPerInt
	jStart := startTwo/basesPerInt
	iEnd := (lenOne+basesPerInt-1)/basesPerInt
	jEnd := (lenTwo+basesPerInt-1)/basesPerInt

	lenOffsetOne := (lenOne % basesPerInt) * bitsPerBase + rightShift
	lenOffsetTwo := (lenTwo % basesPerInt) * bitsPerBase

	seqDiff = seqOne[iStart] >> rightShift
	prevBits = seqOne[iEnd] << leftShift
	seqDiff = seqDiff ^ seqTwo[j]
	seqDiff = seqDiff & (ones >> offsetTwo)

	if iStart+1 >= iEnd && lenOffsetOne != 0 {
		seqDiff = seqDiff | (ones >> lenOffsetOne)
	}
	if jStart+1 >= jEnd && lenOffsetTwo != 0 {
		seqDiff = seqDiff | (ones >> lenOffsetTwo)
	}

	matches = (bits.LeadingZeros64(seqDiff) - offsetTwo) / bitsPerBase

	for i, j = iStart+1, jStart+1; i < (lenOne+basesPerInt-1)/basesPerInt && j < (lenTwo+basesPerInt-1)/basesPerInt; i, j = i+1, j+1 {

		seqDiff = prevBits | (seqOne[i] >> rightShift)

		prevBits = seqOne[i] << leftShift

		seqDiff = seqDiff ^ seqTwo[j]

		if i+1 == (lenOne+basesPerInt-1)/basesPerInt {
			if lenOffsetOne != 0 {
				seqDiff = seqDiff | (ones >> lenOffsetOne)
			}
		}
		if j+1 == (lenTwo+basesPerInt-1)/basesPerInt {
			if lenOffsetTwo != 0 {
				seqDiff = seqDiff | (ones >> lenOffsetTwo)
			}
		}

		matches += bits.LeadingZeros64(seqDiff) / int(bitsPerBase)
	}

	if
	return matches
}*/
