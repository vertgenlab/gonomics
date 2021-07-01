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
		actual = IsSelfOverlapping(b)
		if actual != v.expected {
			t.Errorf("Error in IsSelfOverlapping. Expected: %t. Actual: %t.", v.expected, actual)
		}
	}
}
