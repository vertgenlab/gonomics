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

func TestMaxInt64(t *testing.T) {
	for _, test := range maxTests {
		actual := MaxInt64(int64(test.a), int64(test.b))
		if actual != int64(test.expt) {
			t.Errorf("MaxInt64(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

func TestMaxInt32(t *testing.T) {
	for _, test := range maxTests {
		actual := MaxInt32(int32(test.a), int32(test.b))
		if actual != int32(test.expt) {
			t.Errorf("MaxInt32(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

//like maxTests but for uint types.
var uMaxTests = []struct {
	a uint // first input
	b	uint // second input
	expt uint // expected result
}{
	{1, 2, 2},
	{3, 7, 7},
	{100, 99, 100},
}

func TestMaxUint32(t *testing.T) {
	for _, test := range uMaxTests {
		actual := MaxUint32(uint32(test.a), uint32(test.b))
		if actual != uint32(test.expt) {
			t.Errorf("MaxUint32(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
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

func TestMinInt32(t *testing.T) {
	for _, test := range minTests {
		actual := MinInt32(int32(test.a), int32(test.b))
		if actual != int32(test.expt) {
			t.Errorf("MinInt32(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

func TestMinInt64(t *testing.T) {
	for _, test := range minTests {
		actual := MinInt64(int64(test.a), int64(test.b))
		if actual != int64(test.expt) {
			t.Errorf("MinInt64(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var uMinTests = []struct {
	a    uint // first input
	b    uint // second input
	expt uint // expected result
}{
	{1, 2, 1},
	{3, 7, 3},
	{100, 99, 99},
}

func TestMinUint32(t *testing.T) {
	for _, test := range uMinTests {
		actual := MinUint32(uint32(test.a), uint32(test.b))
		if actual != uint32(test.expt) {
			t.Errorf("MinUint32(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var tripleMinTests = []struct {
	a    int // first input
	b    int // second input
	c 	int //third input
	expt int // expected result
}{
	{1, 2, 10, 1},
	{3, 7, 1,1},
	{100, 99, 2000, 99},
}

func TestTripleMin(t *testing.T) {
	for _, test := range tripleMinTests {
		actual := TripleMin(test.a, test.b, test.c)
		if actual != test.expt {
			t.Errorf("TripleMin(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var tripleMaxTests = []struct {
	a    int // first input
	b    int // second input
	c 	int //third input
	expt int // expected result
}{
	{1, 2, 10, 10},
	{3, 7, 1,7},
	{100, 99, 2000, 2000},
}

func TestTripleMax(t *testing.T) {
	for _, test := range tripleMaxTests {
		actual := TripleMax(test.a, test.b, test.c)
		if actual != test.expt {
			t.Errorf("TripleMax(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var minIntSliceTests = []struct {
	a []int
	expt int
}{
	{[]int{20, 30, 40}, 20},
	{[]int{10, 20, -10, -20, 40}, -20},
}

func TestMinIntSlice(t *testing.T) {
	for _, test := range minIntSliceTests {
		actual := MinIntSlice(test.a)
		if actual != test.expt {
			t.Errorf("MinIntSlice(%v): expected: %d, actual %d.", test.a, test.expt, actual)
		}
	}
}

var maxIntSliceTests = []struct {
	a []int
	expt int
}{
	{[]int{20, 30, 40}, 40},
	{[]int{10, 20, -10, -20, 40}, 40},
}

func TestMaxIntSlice(t *testing.T) {
	for _, test := range maxIntSliceTests {
		actual := MaxIntSlice(test.a)
		if actual != test.expt {
			t.Errorf("MaxIntSlice(%v): expected: %d, actual %d.", test.a, test.expt, actual)
		}
	}
}
