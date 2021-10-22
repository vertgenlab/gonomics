package bed

import "testing"

var TrimTests = []struct {
	InFile       string
	ExpectedFile string
	TrimLeft     int
	TrimRight    int
}{
	{"testdata/testTrim.bed", "testdata/expectedTrim.bed", 10, 10},
}

func TestTrim(t *testing.T) {
	var b, expected []Bed
	for _, v := range TrimTests {
		b = Read(v.InFile)
		expected = Read(v.ExpectedFile)
		Trim(b, v.TrimLeft, v.TrimRight)
		if !AllAreEqual(b, expected) {
			t.Errorf("Error in Trim.")
		}
	}
}
