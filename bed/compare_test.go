package bed

import (
	"testing"
)

var CompareDistanceTests = []struct {
	A 					Bed
	B					Bed
	expectedDistance	int
}{
	//{A: Bed{Chrom:"chr1", ChromStart:10, ChromEnd:20, Name:""}, B: Bed{Chrom:"chr2", ChromStart:30, ChromEnd:40, Name:""}, expectedDistance: -1}, //different chr, only returning panic, DO I REMOVE THIS TEST?
	{A: Bed{Chrom:"chr1", ChromStart:10, ChromEnd:20, Name:""}, B: Bed{Chrom:"chr1", ChromStart:15, ChromEnd:30, Name:""}, expectedDistance: 0}, //test layered coordinates
	{A: Bed{Chrom:"chr1", ChromStart:10, ChromEnd:20, Name:""}, B: Bed{Chrom:"chr1", ChromStart:20, ChromEnd:30, Name:""}, expectedDistance: 1}, //test layered coordinates
	{A: Bed{Chrom:"chr1", ChromStart:10, ChromEnd:20, Name:""}, B: Bed{Chrom:"chr1", ChromStart:30, ChromEnd:40, Name:""}, expectedDistance: 11}, //A upstream
	{A: Bed{Chrom:"chr1", ChromStart:30, ChromEnd:40, Name:""}, B: Bed{Chrom:"chr1", ChromStart:10, ChromEnd:20, Name:""}, expectedDistance: 11}, //B upstream

}

func TestCompareDistance(t *testing.T) {
	var distance int
	for _,v := range CompareDistanceTests {
		distance = CompareDistance(v.A, v.B)
		if distance != v.expectedDistance {
			t.Errorf("Error in CompareDistance. Expected: %v. Actual %v.", v.expectedDistance, distance)
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
