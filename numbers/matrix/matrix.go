package matrix

import (
	"github.com/vertgenlab/gonomics/exception"
	"gonum.org/v1/gonum/mat"
	"log"
	"math"
)

// FractionalSymmetricMatrixExponentiation calculates A^t, where A is a symmetric square matrix *mat.Dense
// and t is a float64 power. Note that t does not have to be an integer, but can be a fractional power.
// A^t is calculated with the identity A^t = exp(t*log(A))
func FractionalSymmetricMatrixExponentiation(inMat *mat.Dense, power float64) *mat.Dense {
	var logInMat *mat.Dense
	var scaledInMat, answer mat.Dense
	logInMat = DenseLogSymmetric(inMat)
	scaledInMat.Scale(power, logInMat)
	answer.Exp(&scaledInMat)
	return &answer
}

// DenseLogSymmetric calculates the logarithm of an input square matrix *mat.Dense using the following identity:
// For a square matrix A, log(A) = V*log(D')*V^{-1}, where A=VDV^-1 is the eigen decomposition of A,
// and D' is the element-wise logarithm of the eigenvalue diagonal matrix.
// Since we take in a symmetric matrix, the logarithm is guaranteed to be real-valued, so we discard the
// imaginary portion of the eigen decomposition.
func DenseLogSymmetric(inMatrix *mat.Dense) *mat.Dense {
	var err error
	var outMatrix mat.Dense
	var complexVectors mat.CDense
	if !IsSymmetric(inMatrix) {
		log.Fatalf("Error: DenseLogSymmetric supports only symmetric input matrices. Input: %v.\n", inMatrix)
	}
	// eigen decomposition
	var eig mat.Eigen
	eigenPassed := eig.Factorize(inMatrix, mat.EigenBoth)
	if !eigenPassed {
		log.Println(mat.Formatted(inMatrix))
		log.Fatalf("Error in eigen decomposition of the above matrix.\n")
	}
	eigenValues := complexSliceToRealSlice(eig.Values(nil))
	for i := range eigenValues {
		eigenValues[i] = math.Log(eigenValues[i])
	}
	eig.VectorsTo(&complexVectors)
	realVectors := complexMatrixToRealMatrix(&complexVectors)
	diagonal := diagonalDenseFromVector(eigenValues)
	outMatrix.Mul(realVectors, diagonal)
	var inverse mat.Dense
	err = inverse.Inverse(realVectors)
	exception.PanicOnErr(err)
	outMatrix.Mul(&outMatrix, &inverse)
	return &outMatrix
}

// complexMatrixToRealMatrix converts a complex-valued matrix *mat.CDense into
// a real-valued matrix *mat.Dense, keeping only the real part and discarding the
// imaginary parts of each value.
func complexMatrixToRealMatrix(complexMat *mat.CDense) *mat.Dense {
	var currRow, currCol int
	rows, cols := complexMat.Dims()
	answer := mat.NewDense(rows, cols, nil)

	for currRow = 0; currRow < rows; currRow++ {
		for currCol = 0; currCol < cols; currCol++ {
			answer.Set(currRow, currCol, real(complexMat.At(currRow, currCol)))
		}
	}
	return answer
}

// diagonalDenseFromVector creates a diagonal dense matrix where the diagonal entries
// are equal to the entries in an input []float64.
func diagonalDenseFromVector(values []float64) *mat.Dense {
	var answer = mat.NewDense(len(values), len(values), nil)
	for currIndex := range values {
		answer.Set(currIndex, currIndex, values[currIndex])
	}
	return answer
}

// IsSymmetric returns true if an input mat.Dense m is symmetric, false otherwise.
// The check is that a symmetric matrix is equal to its transpose.
func IsSymmetric(m *mat.Dense) bool {
	rows, cols := m.Dims() // get the dimensions of the matrix
	if rows != cols {
		return false //only square matrices are symmetric
	}
	transposed := m.T() // transpose the matrix
	return mat.Equal(m, transposed)
}

// complexSliceToRealSlice converts a slice of complex numbers into a slice of floats, taking
// only the real part of the number.
func complexSliceToRealSlice(inSlice []complex128) []float64 {
	var answer = make([]float64, len(inSlice))
	for i := range inSlice {
		answer[i] = real(inSlice[i])
	}
	return answer
}

// Rref returns an input matrix in row reduced echelon form using Gaussian elimination.
func Rref(m [][]float64) [][]float64 {
	//first we make a copy of the input matrix to serve as our answer. This preserves the input matrix in its original form.
	mCopy := make([][]float64, len(m))
	for i := range m {
		mCopy[i] = make([]float64, len(m[i]))
		copy(mCopy[i], m[i])
	}

	var subtractFactor, entry float64
	var leadingCoefficient int = 0
	var i, j int
	for row := 0; row < len(mCopy); row++ {
		if leadingCoefficient >= len(mCopy[0]) {
			return mCopy
		}
		i = row
		for mCopy[i][leadingCoefficient] == 0 {
			i++
			if len(mCopy) == i {
				i = row
				leadingCoefficient++
				if leadingCoefficient == len(mCopy[0]) {
					return mCopy
				}
			}
		}
		mCopy[i], mCopy[row] = mCopy[row], mCopy[i]
		multFactor := 1 / mCopy[row][leadingCoefficient]
		for j = range mCopy[row] {
			mCopy[row][j] *= multFactor
		}
		for i = 0; i < len(mCopy); i++ {
			if i != row {
				subtractFactor = mCopy[i][leadingCoefficient]
				for j, entry = range mCopy[row] {
					mCopy[i][j] -= entry * subtractFactor
				}
			}
		}
		leadingCoefficient++
	}
	return mCopy
}

// approxEqual returns true if the entries of 2 input matrices are approximately equal,
// within some bound of relative precision.
func approxEqual(m1 *mat.Dense, m2 *mat.Dense, precision float64) bool {
	var currRow, currCol int
	rows1, cols1 := m1.Dims()
	rows2, cols2 := m2.Dims()
	if rows1 != rows2 {
		return false
	}
	if cols1 != cols2 {
		return false
	}

	for currRow = 0; currRow < rows1; currRow++ {
		for currCol = 0; currCol < cols1; currCol++ {
			if math.Abs(m1.At(currRow, currCol)-m2.At(currRow, currCol))/math.Max(m1.At(currRow, currCol), m2.At(currRow, currCol)) > precision {
				return false
			}
		}
	}

	return true
}
