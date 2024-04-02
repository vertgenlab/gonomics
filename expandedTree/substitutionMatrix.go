package expandedTree

import (
	"github.com/vertgenlab/gonomics/numbers/matrix"
	"gonum.org/v1/gonum/mat"
)

// PopulateSubstitutionMatrices initializes the SubstitutionMatrix field of every node in an input
// newick tree if the provided node is the root node.
// The substitution matrix for each branch is based on a known unitMatrix [][]float64, which is the
// substitution matrix S measured (often empirically) over a unitBranchLength float64. Thus, for a
// branch of length t, the substitution matrix is S^{t / unitBranchLength}
func PopulateSubstitutionMatrices(node *ETree, unitMatrix [][]float64, unitBranchLength float64) {
	var unitDense, answerDense *mat.Dense
	unitDense = floatMatrixToDenseMatrix(unitMatrix)
	answerDense = matrix.FractionalSymmetricMatrixExponentiation(unitDense, node.BranchLength/unitBranchLength)
	node.SubstitutionMatrix = denseMatrixToFloatMatrix(answerDense)
	if node.Left != nil {
		PopulateSubstitutionMatrices(node.Left, unitMatrix, unitBranchLength)
	}
	if node.Right != nil {
		PopulateSubstitutionMatrices(node.Right, unitMatrix, unitBranchLength)
	}
}

// floatMatrixToDenseMatrix converts an input [][]float64 into a gonum *mat.Dense
func floatMatrixToDenseMatrix(inMat [][]float64) *mat.Dense {
	rows := len(inMat)
	if rows < 1 {
		return mat.NewDense(0, 0, nil)
	}
	cols := len(inMat[0])
	dense := mat.NewDense(rows, cols, nil)
	for currRow := 0; currRow < rows; currRow++ {
		for currCol := 0; currCol < cols; currCol++ {
			dense.Set(currRow, currCol, inMat[currRow][currCol])
		}
	}
	return dense
}

// denseMatrixToFloatMatrix converts a gonum *mat.Dense into a [][]float64
func denseMatrixToFloatMatrix(dense *mat.Dense) [][]float64 {
	rows, cols := dense.Dims()
	var currRow, currCol int

	// Create a new float matrix with the same dimensions as the dense matrix
	floatMat := make([][]float64, rows)
	for i := range floatMat {
		floatMat[i] = make([]float64, cols)
	}

	// Populate the float matrix with values from the dense matrix
	for currRow = 0; currRow < rows; currRow++ {
		for currCol = 0; currCol < cols; currCol++ {
			floatMat[currRow][currCol] = dense.At(currRow, currCol)
		}
	}
	return floatMat
}
