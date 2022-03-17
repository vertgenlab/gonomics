package lift

import (
	"github.com/vertgenlab/gonomics/bed"
	"testing"
)

var OverlapTests = []struct {
	A              bed.Bed
	B              bed.Bed
	expected       bool
	expectedLength int
}{
	{A: bed.Bed{"chr4", 1, 10, "", 0, bed.Positive, 3, nil}, B: bed.Bed{"chr4", 4, 12, "", 0, bed.Positive, 3, nil}, expected: true, expectedLength: 6},
	{A: bed.Bed{"chr5", 1, 10, "", 0, bed.Positive, 3, nil}, B: bed.Bed{"chr4", 4, 12, "", 0, bed.Positive, 3, nil}, expected: false, expectedLength: 0},
	{A: bed.Bed{"chr4", 1, 10, "", 0, bed.Positive, 3, nil}, B: bed.Bed{"chr4", 13, 15, "", 0, bed.Positive, 3, nil}, expected: false, expectedLength: 0},
	{A: bed.Bed{"chr4", 1, 10, "", 0, bed.Positive, 3, nil}, B: bed.Bed{"chr4", 10, 12, "", 0, bed.Positive, 3, nil}, expected: false, expectedLength: 0},
}

func TestOverlap(t *testing.T) {
	var actual bool
	for _, v := range OverlapTests {
		actual = overlap(v.A, v.B)
		if actual != v.expected {
			t.Errorf("Error in Overlap. Expected: %t. Actual: %t.", v.expected, actual)
		}
	}
}

var OverlapCountTests = []struct {
	elements1File string
	elements2File string
	expected      int
	expectedSum   int
}{
	{"testdata/elements1.bed", "testdata/elements2.bed", 1, 2},
}

func TestOverlapCount(t *testing.T) {
	var elements1, elements2 []Lift
	var actual int

	for _, v := range OverlapCountTests {
		elements1 = GoRead(v.elements1File)
		elements2 = GoRead(v.elements2File)
		actual = OverlapCount(elements1, elements2)
		if actual != v.expected {
			t.Errorf("Error in OverlapCount. Expected: %d. Actual: %d.", v.expected, actual)
		}
	}
}

func TestOverlapLength(t *testing.T) {
	var actual int
	for _, v := range OverlapTests {
		actual = overlapLength(v.A, v.B)
		if actual != v.expectedLength {
			t.Errorf("Error in OverlapLength. Expected: %d. Actual: %d.", v.expectedLength, actual)
		}
	}
}

func TestOverlapLengthSum(t *testing.T) {
	var elements1, elements2 []Lift
	var actual int
	for _, v := range OverlapCountTests {
		elements1 = GoRead(v.elements1File)
		elements2 = GoRead(v.elements2File)
		actual = overlapLengthSum(elements1, elements2)
		if actual != v.expectedSum {
			t.Errorf("Error in OverlapLengthSum. Expected: %d. Actual: %d.", v.expectedSum, actual)
		}
	}
}

var SortTests = []struct {
	inputFile                     string
	expectedByCoordFile           string
	expectedBySizeFile            string
	expectedByChromEndByChromFile string
}{
	{"testdata/sortInput.bed", "testdata/expectedSortByCoord.bed", "testdata/expectedSortBySize.bed", "testdata/expectedSortByChromEndByChrom.bed"},
}

func TestSortByCoord(t *testing.T) {
	var input, expectedCoord []Lift
	for _, v := range SortTests {
		input = GoRead(v.inputFile)
		expectedCoord = GoRead(v.expectedByCoordFile)
		SortByCoord(input)
		if !allAreEqual(input, expectedCoord) {
			t.Errorf("Error in SortByCoord.")
		}
	}
}
