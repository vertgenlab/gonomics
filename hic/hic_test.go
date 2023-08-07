package hic

import (
	"testing"
)

var straw1 = Straw{Bin1Start: 5000, Bin2Start: 1000, ContactScore: 5}
var straw2 = Straw{Bin1Start: 2000, Bin2Start: 10000, ContactScore: 10}
var straws = []Straw{straw1, straw2}

var readStrawTests = []struct {
	filename string
	data     []Straw
}{
	{"testdata/strawTestFile.straw", straws},
}

func TestRead(t *testing.T) {
	var actual []Straw
	var ch <-chan Straw
	var index int
	for _, test := range readStrawTests {
		actual = Read(test.filename)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The test data does not match the data in the file %s", test.filename)
		}

		index = 0
		ch = GoReadToChan(test.filename)
		for i := range ch {
			if !Equal(i, test.data[index]) {
				t.Errorf("Error: GoReadToChan output does not match expected.\n")
			}
			index++
		}

	}
}
