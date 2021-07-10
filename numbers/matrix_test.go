package numbers

import (
	"testing"
	"fmt"
)

var RrefTests = []struct {
	input [][]float64
	expected [][]float64
}{
	{[][]float64{{1, 1, 7},{1, 2, 11}}, [][]float64{{1, 0, 3},{0, 1, 4}}},
	{[][]float64{{1, 2, -1, -4},{2, 3, -1, -11}, {-2, 0, -3, 22}}, [][]float64{{1, 0, 0, -8},{0, 1, 0, 1}, {0, 0, 1, -2}}},
}

func TestRref(t *testing.T) {
	var output [][]float64
	for _, v := range RrefTests {
		output = Rref(v.input)
		if !equalMatrix(output, v.expected) {
			fmt.Println(output)
			t.Errorf("Error in Rref.")
		}
	}
}

func equalMatrix(a [][]float64, b [][]float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if len(a[i]) != len(b[i]) {
			return false
		}
		for j := range a[i] {
			if a[i][j] != b[i][j] {
				return false
			}
		}
	}
	return true
}