package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

func TestView(t *testing.T) {
	// Test Data
	alpha := dna.StringToBases("TTTTTTTTTTTTTTTTAGC")
	beta := dna.StringToBases("ATTTTTTTAGC")
	operations := []cigar.Cigar{
		{RunLength: 10, Op: '='},
	}
	startI, endI, startJ, endJ := 9, 19, 1, 11

	// Expected Output
	expectedAlignment := "TTTTTTTTT-TTTTTTTAGC\n---------ATTTTTTTAGC\n"

	// Function Call
	alignment := View(alpha, beta, operations, startI, endI, startJ, endJ)

	t.Logf("\n\n%s\n", alignment)

	// Assertion
	if alignment != expectedAlignment {
		t.Errorf("Alignment mismatch:\nExpected:\n%s\nGot:\n%s", expectedAlignment, alignment)
	}
}
