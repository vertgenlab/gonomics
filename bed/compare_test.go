package bed

import (
	"testing"
)

var CompareDistanceTests = []struct {
	A                Bed
	B                Bed
	expectedDistance int
	expectedError    bool
}{
	{A: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: ""}, B: Bed{Chrom: "chr2", ChromStart: 21, ChromEnd: 30, Name: ""}, expectedDistance: -1, expectedError: true},  //test if input beds are on diff chrom
	{A: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: ""}, B: Bed{Chrom: "chr1", ChromStart: 15, ChromEnd: 30, Name: ""}, expectedDistance: 0, expectedError: false},  //test layered coordinates
	{A: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: ""}, B: Bed{Chrom: "chr1", ChromStart: 20, ChromEnd: 30, Name: ""}, expectedDistance: 1, expectedError: false},  //test layered coordinates
	{A: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: ""}, B: Bed{Chrom: "chr1", ChromStart: 30, ChromEnd: 40, Name: ""}, expectedDistance: 11, expectedError: false}, //A upstream
	{A: Bed{Chrom: "chr1", ChromStart: 30, ChromEnd: 40, Name: ""}, B: Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Name: ""}, expectedDistance: 11, expectedError: false}, //B upstream
}

func TestCompareDistance(t *testing.T) {
	var distance int
	var err error
	var boolForError bool //to help test if the error returned is as expected
	for _, v := range CompareDistanceTests {
		distance, err = MinimumDistance(v.A, v.B)
		if distance != v.expectedDistance {
			t.Errorf("Error in CompareDistance. Expected: %v. Actual %v.", v.expectedDistance, distance)
		}
		if err == nil {
			boolForError = false
		} else {
			boolForError = true
		}
		if boolForError != v.expectedError {
			t.Errorf("Error in CompareDistance. Expected error: %t. Actual error: %t", v.expectedError, boolForError)
		}
	}
}

var OverlapTests = []struct {
	A              Bed
	B              Bed
	expected       bool
	expectedLength int
}{
	{A: Bed{"chr4", 1, 10, "", 0, Positive, 3, nil}, B: Bed{"chr4", 4, 12, "", 0, Positive, 3, nil}, expected: true, expectedLength: 6},
	{A: Bed{"chr5", 1, 10, "", 0, Positive, 3, nil}, B: Bed{"chr4", 4, 12, "", 0, Positive, 3, nil}, expected: false, expectedLength: 0},
	{A: Bed{"chr4", 1, 10, "", 0, Positive, 3, nil}, B: Bed{"chr4", 13, 15, "", 0, Positive, 3, nil}, expected: false, expectedLength: 0},
	{A: Bed{"chr4", 1, 10, "", 0, Positive, 3, nil}, B: Bed{"chr4", 10, 12, "", 0, Positive, 3, nil}, expected: false, expectedLength: 0},
}

func TestOverlap(t *testing.T) {
	var actual bool
	for _, v := range OverlapTests {
		actual = Overlap(v.A, v.B)
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
	var elements1, elements2 []Bed
	var actual int

	for _, v := range OverlapCountTests {
		elements1 = Read(v.elements1File)
		elements2 = Read(v.elements2File)
		actual = OverlapCount(elements1, elements2)
		if actual != v.expected {
			t.Errorf("Error in OverlapCount. Expected: %d. Actual: %d.", v.expected, actual)
		}
	}
}

func TestOverlapLength(t *testing.T) {
	var actual int
	for _, v := range OverlapTests {
		actual = OverlapLength(v.A, v.B)
		if actual != v.expectedLength {
			t.Errorf("Error in OverlapLength. Expected: %d. Actual: %d.", v.expectedLength, actual)
		}
	}
}

func TestOverlapLengthSum(t *testing.T) {
	var elements1, elements2 []Bed
	var actual int
	for _, v := range OverlapCountTests {
		elements1 = Read(v.elements1File)
		elements2 = Read(v.elements2File)
		actual = OverlapLengthSum(elements1, elements2)
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
	var input, expectedCoord []Bed
	for _, v := range SortTests {
		input = Read(v.inputFile)
		expectedCoord = Read(v.expectedByCoordFile)
		SortByCoord(input)
		if !AllAreEqual(input, expectedCoord) {
			t.Errorf("Error in SortByCoord.")
		}
	}
}

func TestSortBySize(t *testing.T) {
	var input, expectedSize []Bed
	for _, v := range SortTests {
		input = Read(v.inputFile)
		expectedSize = Read(v.expectedBySizeFile)
		SortBySize(input)
		if !AllAreEqual(input, expectedSize) {
			t.Errorf("Error in SortBySize.")
		}
	}
}

func TestSortByChromEndByChrom(t *testing.T) {
	var input, expectedChromEnd []Bed
	for _, v := range SortTests {
		input = Read(v.inputFile)
		expectedChromEnd = Read(v.expectedByChromEndByChromFile)
		SortByChromEndByChrom(input)
		if !AllAreEqual(input, expectedChromEnd) {
			t.Errorf("Error in SortByChromEndbyChrom.")
		}
	}
}

var Equal1 Bed = Bed{Chrom: "chr1", ChromStart: 100, ChromEnd: 200, Name: "First", Score: 1, Strand: '+', FieldsInitialized: 6}
var Equal2 Bed = Bed{Chrom: "chr2", ChromStart: 400, ChromEnd: 900, Name: "Second", Score: 5, Strand: '-', FieldsInitialized: 6}

var EqualTests = []struct {
	a        Bed
	b        Bed
	Expected bool
}{
	{a: Equal1, b: Equal2, Expected: false},
	{a: Equal1, b: Equal1, Expected: true},
}

func TestEqual(t *testing.T) {
	var actual bool
	for _, v := range EqualTests {
		actual = Equal(v.a, v.b)
		if actual != v.Expected {
			t.Errorf("Error in bed package Equal function.")
		}
	}
}

var AllEqualTests = []struct {
	a        []Bed
	b        []Bed
	Expected bool
}{
	{a: []Bed{Equal1, Equal2, Equal2}, b: []Bed{Equal1}, Expected: false},
	{a: []Bed{Equal2, Equal1}, b: []Bed{Equal2, Equal1}, Expected: true},
	{a: []Bed{Equal2, Equal2, Equal1, Equal1}, b: []Bed{Equal2, Equal1, Equal1, Equal1}, Expected: false},
}

func TestAllAreEqual(t *testing.T) {
	var actual bool
	for _, v := range AllEqualTests {
		actual = AllAreEqual(v.a, v.b)
		if actual != v.Expected {
			t.Errorf("Error in bed package AllAreEqual function.")
		}
	}
}
