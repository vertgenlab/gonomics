package dna

import (
	"math"
	"testing"
)

func TestMeltingTemp(t *testing.T) {
	var testSeqs [][]Base = [][]Base{StringToBases("ATCGTGACTGA"), StringToBases("GTCGTGATTCTGC"), StringToBases("GTCGTTAGATTCTGT"), StringToBases("GCTGCGAATTCGCAGC")}
	var expectedTm []float64 = []float64{32.4608090067, 41.6641715041, 41.0485726487, 55.4258364707}

	for i := range testSeqs {
		if math.Abs(expectedTm[i]-MeltingTemp(testSeqs[i])) > 1e-6 {
			t.Errorf("Error when calculating the Tm for: %s. Output: %.10f  Expected: %.10f", BasesToString(testSeqs[i]), MeltingTemp(testSeqs[i]), expectedTm[i])
		}
	}
}
