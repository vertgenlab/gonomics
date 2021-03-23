package numbers

import (
	"testing"
)

//TODO: Testing functions for min/max functions for other data types.

var maxTests = []struct {
	a    int // first input
	b    int // second input
	expt int // expected result
}{
	{1, 2, 2},
	{-1, 2, 2},
	{-50, 2, 2},
	{3, 7, 7},
	{3, -7, 3},
	{-3, -7, -3},
	{100, 99, 100},
	{100, -101, 100},
}

func TestMax(t *testing.T) {
	for _, test := range maxTests {
		actual := Max(test.a, test.b)
		if actual != test.expt {
			t.Errorf("Max(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var minTests = []struct {
	a    int // first input
	b    int // second input
	expt int // expected result
}{
	{1, 2, 1},
	{-1, 2, -1},
	{-50, 2, -50},
	{3, 7, 3},
	{3, -7, -7},
	{-3, -7, -7},
	{100, 99, 99},
	{100, -101, -101},
}

func TestMin(t *testing.T) {
	for _, test := range minTests {
		actual := Min(test.a, test.b)
		if actual != test.expt {
			t.Errorf("Min(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}
