package cigar

import "testing"

func TestIsByteCigarEqual(t *testing.T) {
	// Test case 1: Equal ByteCigars
	a := ByteCigar{RunLen: 3, Op: 'M'}
	b := ByteCigar{RunLen: 3, Op: 'M'}
	if !isByteCigarEqual(a, b) {
		t.Errorf("Expected ByteCigars to be equal, but they are not")
	}

	// Test case 2: Different RunLen
	c := ByteCigar{RunLen: 5, Op: 'M'}
	if isByteCigarEqual(a, c) {
		t.Errorf("Expected ByteCigars to be different, but they are equal")
	}

	// Test case 3: Different Op
	d := ByteCigar{RunLen: 3, Op: 'I'}
	if isByteCigarEqual(a, d) {
		t.Errorf("Expected ByteCigars to be different, but they are equal")
	}

	// Test slice

	alpha := []ByteCigar{{RunLen: 1, Op: 'M'}, {RunLen: 1, Op: 'D'}, {RunLen: 1, Op: 'M'}}
	beta := []ByteCigar{{RunLen: 1, Op: 'M'}, {RunLen: 1, Op: 'D'}, {RunLen: 1, Op: 'M'}}

	if !EqualByteCigar(alpha, beta) {
		t.Errorf("Expected ByteCigars to be different, but they are equal")
	}
}
