package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

func TestGetTargetBases(t *testing.T) {
	// Test case 1: direction is left
	n := &Node{
		Seq: dna.StringToBases("ACGT"),
	}
	extension := len(n.Seq)
	position := len(n.Seq) - 1
	direction := left

	expected := "ACG"
	result := dna.BasesToString(getTargetBases(n, extension, position, []dna.Base{}, []dna.Base{}, direction))

	if result != expected {
		t.Errorf("Test case 1 failed. Expected: %v, got: %v", expected, result)
	}

	// Test case 2: direction is right
	n = &Node{
		Seq: dna.StringToBases("ACGT"),
	}
	extension = 5
	position = 0

	direction = right

	expected = "ACGT"
	result = dna.BasesToString(getTargetBases(n, extension, position, []dna.Base{}, []dna.Base{}, direction))

	if result != expected {
		t.Errorf("Test case 2 failed. Expected: %v, got: %v", expected, result)
	}

	// Add more test cases here...
}
