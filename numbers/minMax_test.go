package numbers

import (
	"testing"
)

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

var tripleMinTests = []struct {
	a    int // first input
	b    int // second input
	c    int //third input
	expt int // expected result
}{
	{1, 2, 10, 1},
	{3, 7, 1, 1},
	{100, 99, 2000, 99},
}

func TestTripleMin(t *testing.T) {
	for _, test := range tripleMinTests {
		actual := MinMany(test.a, test.b, test.c)
		if actual != test.expt {
			t.Errorf("TripleMin(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var tripleMaxTests = []struct {
	a    int // first input
	b    int // second input
	c    int //third input
	expt int // expected result
}{
	{1, 2, 10, 10},
	{3, 7, 1, 7},
	{100, 99, 2000, 2000},
}

func TestTripleMax(t *testing.T) {
	for _, test := range tripleMaxTests {
		actual := MaxMany(test.a, test.b, test.c)
		if actual != test.expt {
			t.Errorf("TripleMax(%d, %d): expected %d, actual %d", test.a, test.b, test.expt, actual)
		}
	}
}

var minIntSliceTests = []struct {
	a    []int
	expt int
}{
	{[]int{20, 30, 40}, 20},
	{[]int{10, 20, -10, -20, 40}, -20},
}

func TestMinIntSlice(t *testing.T) {
	for _, test := range minIntSliceTests {
		actual := MinMany(test.a...)
		if actual != test.expt {
			t.Errorf("MinIntSlice(%v): expected: %d, actual %d.", test.a, test.expt, actual)
		}
	}
}

var maxIntSliceTests = []struct {
	a    []int
	expt int
}{
	{[]int{20, 30, 40}, 40},
	{[]int{10, 20, -10, -20, 40}, 40},
}

func TestMaxIntSlice(t *testing.T) {
	for _, test := range maxIntSliceTests {
		actual := MaxMany(test.a...)
		if actual != test.expt {
			t.Errorf("MaxIntSlice(%v): expected: %d, actual %d.", test.a, test.expt, actual)
		}
	}
}

func BenchmarkTwoMin(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Min(2, 1)
	}
}

func BenchmarkAnyMin(b *testing.B) {
	for i := 0; i < b.N; i++ {
		MinMany(2, 1)
	}
}
