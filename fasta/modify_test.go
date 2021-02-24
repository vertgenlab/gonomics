package fasta

import "testing"

func TestRemove(t *testing.T) {
	for _, test := range readWriteTests {
		rmCopy := *test.data[1]
		expected := []*Fasta{test.data[0], test.data[2]}
		actual := Remove(test.data, 1)
		if !AllAreEqual(expected, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		// add back so other tests are not broken
		test.data[1] = &rmCopy
	}
}

func TestReverseComplement(t *testing.T) {
	for _, test := range allRevCompTests {
		ReverseComplementAll(test.input)
		if !AllAreEqual(test.input, test.expected) {
			t.Errorf("Expected reverse complement to give %v, but got %v.", test.input, test.expected)
		}
	}
}
