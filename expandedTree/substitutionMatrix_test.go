package expandedTree

import (
	"fmt"
	"math"
	"testing"
)

func initializeMyRoot() *ETree {
	var myRoot, firstLeft, firstRight, secondLeft, secondRight ETree
	myRoot = ETree{
		Name:  "myRoot",
		Up:    nil,
		Left:  &firstLeft,
		Right: &firstRight,
	}
	firstLeft = ETree{
		Name:         "firstLeft",
		BranchLength: 0.7,
		Up:           &myRoot,
		Left:         &secondLeft,
		Right:        &secondRight,
	}
	firstRight = ETree{
		Name:         "firstRight",
		BranchLength: 1.4,
		Up:           &myRoot,
		Left:         nil,
		Right:        nil,
	}
	secondLeft = ETree{
		Name:         "secondLeft",
		BranchLength: 0,
		Up:           &firstLeft,
		Left:         nil,
		Right:        nil,
	}
	secondRight = ETree{
		Name:         "secondRight",
		BranchLength: 500,
		Up:           &firstLeft,
		Left:         nil,
		Right:        nil,
	}
	return &myRoot
}

var PopulateSubstitutionMatricesTests = []struct {
	Root             *ETree
	UnitMatrix       [][]float64
	UnitBranchLength float64
	Precision        float64
	ExpectedMatrices map[string][][]float64
}{
	{Root: initializeMyRoot(),
		UnitMatrix: [][]float64{
			{0.91, 0.03, 0.03, 0.03},
			{0.03, 0.91, 0.03, 0.03},
			{0.03, 0.03, 0.91, 0.03},
			{0.03, 0.03, 0.03, 0.91},
		},
		UnitBranchLength: 1,
		Precision:        1e-3,
		ExpectedMatrices: map[string][][]float64{
			"firstLeft": {
				{0.936, 0.0214, 0.0214, 0.0214},
				{0.0214, 0.936, 0.0214, 0.0214},
				{0.0214, 0.0214, 0.936, 0.0214},
				{0.0214, 0.0214, 0.0214, 0.936},
			},
			"firstRight": {
				{0.877, 0.041, 0.041, 0.041},
				{0.041, 0.877, 0.041, 0.041},
				{0.041, 0.041, 0.877, 0.041},
				{0.041, 0.041, 0.041, 0.877},
			},
			"secondLeft": {
				{1, 0, 0, 0},
				{0, 1, 0, 0},
				{0, 0, 1, 0},
				{0, 0, 0, 1},
			},
			"secondRight": {
				{0.25, 0.25, 0.25, 0.25},
				{0.25, 0.25, 0.25, 0.25},
				{0.25, 0.25, 0.25, 0.25},
				{0.25, 0.25, 0.25, 0.25},
			},
		},
	},
}

func TestPopulateSubstitutionMatrices(t *testing.T) {
	var nodes []*ETree
	var currIndex int
	for _, v := range PopulateSubstitutionMatricesTests {
		PopulateSubstitutionMatrices(v.Root, v.UnitMatrix, v.UnitBranchLength)
		nodes = GetTree(v.Root)
		for currIndex = range nodes {
			if nodes[currIndex].Up != nil && !approxEqual(nodes[currIndex].SubstitutionMatrix, v.ExpectedMatrices[nodes[currIndex].Name], v.Precision) {
				fmt.Printf("Node ID: %v.\n", nodes[currIndex].Name)
				fmt.Printf("SubstitutionMatrix:\n%v\n", nodes[currIndex].SubstitutionMatrix)
				t.Errorf("Error: PopulateSubstitutionMatrices did not produce the expected result.\n")
			}
		}
	}
}

// approxEqual returns true if the entries of 2 input matrices are approximately equal,
// within some bound of relative precision.
func approxEqual(m1 [][]float64, m2 [][]float64, precision float64) bool {
	var currRow, currCol int
	rows1, cols1 := len(m1), len(m1[0])
	rows2, cols2 := len(m2), len(m2[0])
	if rows1 != rows2 {
		return false
	}
	if cols1 != cols2 {
		return false
	}

	for currRow = 0; currRow < rows1; currRow++ {
		for currCol = 0; currCol < cols1; currCol++ {
			if math.Abs(m1[currRow][currCol]-m2[currRow][currCol])/math.Max(m1[currRow][currCol], m2[currRow][currCol]) > precision {
				return false
			}
		}
	}

	return true
}
