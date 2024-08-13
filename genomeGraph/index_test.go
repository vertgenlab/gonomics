package genomeGraph

import (
	"math/rand"
	"testing"
)


func TestNumberToChromAndPos(t *testing.T) {
    testCases := []struct {
        code       uint64
        expChrom   int64
        expPos     int64
    }{
        {0b0000000000000000000000000000000000000000000000000000000000001010, 0, 10},       // Chrom 0, Pos 10
        {0b000000000000000000000000000000010000000000000000000000001111111, 0, 2147483775},     // Chrom 1, Pos 127 (max for 32-bit pos)
        {0b000000000000000000000000000000100000000000000000000000000000001, 1, 1},       // Chrom 2, Pos 1 (min for 32-bit pos)
        {0b0100000000000000000000000000000000000000000000000000000000000000, 1073741824, 0},    // Chrom 16, Pos 0
    }

    for _, tc := range testCases {
		chrom, pos := numberToChromAndPos(tc.code)

        if chrom != tc.expChrom || pos != tc.expPos {
            t.Errorf("NumberToChromAndPos(%b) = (%d, %d), want (%d, %d)", tc.code, chrom, pos, tc.expChrom, tc.expPos)
        }
		if ChromAndPosToNumber(int(chrom), int(pos)) != tc.code {
			t.Errorf("ChromAndPosToNumber(%d, %d) != %v\n", chrom, pos, tc.code)
		}
    }
}

func BenchmarkOldNumberToChromAndPos(b *testing.B) {
    // Generate random test data
    randCodes := make([]uint64, b.N) // Preallocate to avoid allocations in the loop
    for i := 0; i < b.N; i++ {
        randCodes[i] = rand.Uint64() 
    }

    b.ResetTimer() // Start timing after data generation

    for i := 0; i < b.N; i++ {
        numberToChromAndPos(randCodes[i]) 
    }
	b.ReportAllocs()
}
