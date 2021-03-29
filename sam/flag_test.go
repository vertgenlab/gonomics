package sam

import (
	"fmt"
	"testing"
)

const (
	paired int = 1
	properlyAligned int = 2
	unmapped int = 4
	mateUnmapped int = 8
	posStrand int = 16
	mateIsPosStrand int = 32
	forwardRead int = 64
	reverseRead int = 128
	notPrimeryAlign int = 256
	readFailsQc int = 512
	isDuplicate int = 1024
	suppAlign int = 2048
)

var testMap = map[int]func(sam Aln) bool {
	1: IsPaired,
	2: ProperlyAligned,
	4: IsUnmapped,
	8: MateIsUnmapped,
	16: IsPosStrand,
	32: MateIsPosStrand,
	64: IsForwardRead,
	128: IsReverseRead,
	256: IsNotPrimaryAlign,
	512: ReadFailsQc,
	1024: IsDuplicate,
	2048: IsSupplementaryAlign,
}

func flagPasses(r Aln, trueBits map[int]bool) bool {
	for testedBit, testFunc := range testMap {
		if testFunc(r) != trueBits[testedBit] {
			fmt.Println(r, trueBits)
			fmt.Println(testedBit)
			fmt.Println(IsPaired(r))
			return false
		}
	}
	return true
}

func makeTestRead(trueBits []int) (Aln, map[int]bool) {
	var flag uint16
	truthMap := make(map[int]bool)
	for _, bit := range trueBits {
		flag += uint16(bit)
		truthMap[bit] = true
	}
	return Aln{Flag: flag}, truthMap
}

func testFlags(toTest []int) bool {
	testRead, truthMap := makeTestRead(toTest)
	return flagPasses(testRead, truthMap)
}

func TestFlagFuncs(t *testing.T) {

	if !testFlags([]int{paired, properlyAligned, unmapped, mateUnmapped,
		posStrand, mateIsPosStrand, forwardRead, reverseRead, notPrimeryAlign, readFailsQc,
	    isDuplicate, suppAlign}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{unmapped, mateUnmapped, posStrand, mateIsPosStrand}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{suppAlign}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{reverseRead, notPrimeryAlign, readFailsQc, isDuplicate, suppAlign}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{mateUnmapped, posStrand, forwardRead}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{paired, properlyAligned, forwardRead, reverseRead, readFailsQc, isDuplicate}) {
		t.Errorf("error with flag funcs")
	}

	if !testFlags([]int{paired, unmapped, mateUnmapped, posStrand}) {
		t.Errorf("error with flag funcs")
	}
}
