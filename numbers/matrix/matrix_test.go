package matrix

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
	"testing"
)

var FractionalSymmetricMatrixExponentiationTests = []struct {
	InMat     *mat.Dense
	Power     float64
	Expected  *mat.Dense
	Precision float64
}{
	{InMat: mat.NewDense(4, 4, []float64{
		0.91, 0.03, 0.03, 0.03,
		0.03, 0.91, 0.03, 0.03,
		0.03, 0.03, 0.91, 0.03,
		0.03, 0.03, 0.03, 0.91,
	}),
		Power: 0, //a matrix raised to the zero power is the identity matrix
		Expected: mat.NewDense(4, 4, []float64{
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1,
		}),
		Precision: 1e-6,
	},
	{InMat: mat.NewDense(4, 4, []float64{
		0.91, 0.03, 0.03, 0.03,
		0.03, 0.91, 0.03, 0.03,
		0.03, 0.03, 0.91, 0.03,
		0.03, 0.03, 0.03, 0.91,
	}),
		Power: 500, //substitution matrices approach a flat matrix as time approaches infinity
		Expected: mat.NewDense(4, 4, []float64{
			0.25, 0.25, 0.25, 0.25,
			0.25, 0.25, 0.25, 0.25,
			0.25, 0.25, 0.25, 0.25,
			0.25, 0.25, 0.25, 0.25,
		}),
		Precision: 1e-6,
	},
	{InMat: mat.NewDense(4, 4, []float64{
		0.91, 0.03, 0.03, 0.03,
		0.03, 0.91, 0.03, 0.03,
		0.03, 0.03, 0.91, 0.03,
		0.03, 0.03, 0.03, 0.91,
	}),
		Power: 0.5, //square root
		Expected: mat.NewDense(4, 4, []float64{
			0.954, 0.0154, 0.0154, 0.0154,
			0.0154, 0.954, 0.0154, 0.0154,
			0.0154, 0.0154, 0.954, 0.0154,
			0.0154, 0.0154, 0.0154, 0.954,
		}),
		Precision: 1e-2,
	},
}

func TestFractionalSymmetricMatrixExponentiation(t *testing.T) {
	for _, v := range FractionalSymmetricMatrixExponentiationTests {
		var observed *mat.Dense
		observed = FractionalSymmetricMatrixExponentiation(v.InMat, v.Power)
		if !approxEqual(observed, v.Expected, v.Precision) {
			fmt.Println(mat.Formatted(observed))
			t.Errorf("Error: FractionalSymmetricMatrixExponentiation output was not as expected.")
		}
	}
}

/*
These results can be validated with the independent implementation
from scipy:
import numpy as np
from scipy.linalg import logm
A = np.array([[4, 1], [1, 4]])
log_A = logm(A)
*/
var DenseLogSymmetricTests = []struct {
	Input    *mat.Dense
	Expected *mat.Dense
}{
	{Input: mat.NewDense(2, 2, []float64{
		4, 1,
		1, 4}),
		Expected: mat.NewDense(2, 2, []float64{
			1.3540251005511048, 0.25541281188299536,
			0.25541281188299536, 1.3540251005511048}),
	},
}

func TestDenseLogSymmetric(t *testing.T) {
	for _, v := range DenseLogSymmetricTests {
		var output mat.Dense // we have to allocate in the loop to reset the dimensions
		output = *DenseLogSymmetric(v.Input)
		if !mat.Equal(&output, v.Expected) {
			fmt.Println(mat.Formatted(&output))
			t.Errorf("Error: DenseLogSymmetric output not as expected.")
		}
	}
}

var ComplexMatrixToRealMatrixTests = []struct {
	InMat       *mat.CDense
	ExpectedMat *mat.Dense
}{
	{InMat: mat.NewCDense(4, 4, []complex128{
		4 + 3i, 3, -1 - 2i, 4,
		6 + 2i, 1i, -2, 4,
		5, 5 + 3.4i, 3.5 + 3.2i, 4,
		3, 3, 3, 3i,
	}),
		ExpectedMat: mat.NewDense(4, 4, []float64{
			4, 3, -1, 4,
			6, 0, -2, 4,
			5, 5, 3.5, 4,
			3, 3, 3, 0,
		})},
}

func TestComplexMatrixToRealMatrix(t *testing.T) {
	var observed *mat.Dense
	for _, v := range ComplexMatrixToRealMatrixTests {
		observed = complexMatrixToRealMatrix(v.InMat)
		if !mat.Equal(observed, v.ExpectedMat) {
			t.Errorf("Error: complexMtrixToRealMatrix not as expected.")
		}
	}
}

var DiagonalDenseFromVectorTests = []struct {
	Values   []float64
	Expected *mat.Dense
}{
	{Values: []float64{5, 2, 3, 3},
		Expected: mat.NewDense(4, 4, []float64{
			5, 0, 0, 0,
			0, 2, 0, 0,
			0, 0, 3, 0,
			0, 0, 0, 3})},
}

func TestDiagonalDenseFromVector(t *testing.T) {
	var observed *mat.Dense
	for _, v := range DiagonalDenseFromVectorTests {
		observed = diagonalDenseFromVector(v.Values)
		if !mat.Equal(v.Expected, observed) {
			t.Errorf("Error: DiagonalDenseFromVector not as expected.\n")
		}
	}
}

var IsSymmetricTests = []struct {
	InMat    *mat.Dense
	Expected bool
}{
	{InMat: mat.NewDense(2, 2, []float64{
		4, 1,
		1, 4}),
		Expected: true},
	{InMat: mat.NewDense(3, 3, []float64{
		10, 1, 1,
		1, 10, 1,
		1, 1, 10}),
		Expected: true},
	{InMat: mat.NewDense(3, 3, []float64{
		10, 1, 2,
		1, 10, 1,
		1, 1, 10}),
		Expected: false},
}

func TestIsSymmetric(t *testing.T) {
	var observed bool
	for _, v := range IsSymmetricTests {
		observed = IsSymmetric(v.InMat)
		if observed != v.Expected {
			t.Errorf("Error: IsSymmetric did not have the expected result.")
		}
	}
}

var RrefTests = []struct {
	input    [][]float64
	expected [][]float64
}{
	{[][]float64{{1, 1, 7}, {1, 2, 11}}, [][]float64{{1, 0, 3}, {0, 1, 4}}},
	{[][]float64{{1, 2, -1, -4}, {2, 3, -1, -11}, {-2, 0, -3, 22}}, [][]float64{{1, 0, 0, -8}, {0, 1, 0, 1}, {0, 0, 1, -2}}},
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
