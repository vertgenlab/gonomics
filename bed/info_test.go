package bed

import "testing"

var TotalSizeTests = []struct {
	filename string
	expected int
}{
	{"testdata/elements1.bed", 460},
	{"testdata/elements2.bed", 6},
}

func TestTotalSize(t *testing.T) {
	var b []Bed
	var actual int
	for _, v := range TotalSizeTests {
		b = Read(v.filename)
		actual = TotalSize(b)
		if actual != v.expected {
			t.Errorf("Error in TotalSize. Expected: %d. Actual: %d.", v.expected, actual)
		}
	}
}

var SelfOverlappingTests = []struct {
	filename string
	expected bool
}{
	{"testdata/elements1.bed", false},
	{"testdata/selfOverlap.bed", true},
}

func TestIsSelfOverlapping(t *testing.T) {
	var b []Bed
	var actual bool
	for _, v := range SelfOverlappingTests {
		b = Read(v.filename)
		actual = IsSelfOverlapping(b, 0)
		if actual != v.expected {
			t.Errorf("Error in IsSelfOverlapping. Expected: %t. Actual: %t.", v.expected, actual)
		}
	}
}

func TestOverlapSize(t *testing.T) {
	var answer []int
	var expected []int = []int{3, 2, 0}
	var b Bed = Bed{Chrom: "chr1", ChromStart: 12, ChromEnd: 15}
	var bedA []Bed = []Bed{
		{Chrom: "chr1", ChromStart: 12, ChromEnd: 15},
		{Chrom: "chr1", ChromStart: 13, ChromEnd: 15},
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 500}}
	
	for _, bd := range bedA {
		answer = append(answer, OverlapSize(bd, b))
	}
	for i := range answer {
		if answer[i] == expected[i] {
			continue
		}
		t.Errorf("Error in OverlapSize. Expected: %v. Actual: %v", expected, answer)
	}
}
