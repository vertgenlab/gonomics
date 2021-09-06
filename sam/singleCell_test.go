package sam

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var singleCellTests = []struct {
	InFile      string
	ExpectedBx  []string
	ExpectedUmi []string
}{
	{"testdata/singleCell.sam", []string{"AGTGCACGTGGAGGCT", "AGTGCACGTGGAGGCT", "AGTGCACGTGGAGGCT", "GGACAGTTAGACAAGT", "GCTAGGACAATGAGTA", "CAGGACTGAGTGACCA"}, []string{"AGGTGGACGCTA", "AGGTGGACGCTA", "AGGTGGACGCTA", "GGACAGTGATGA", "CAGTGAGCAGTA", "AGGTGGACGCTA"}},
}

func TestToSingleCellAlignment(t *testing.T) {
	var scAln SingleCellAlignment
	for _, v := range singleCellTests {
		records, _ := Read(v.InFile)
		for i := range records {
			scAln = ToSingleCellAlignment(records[i])
			if dna.BasesToString(scAln.Umi) != v.ExpectedUmi[i] {
				t.Errorf("Error in ToSingleCellAlignment. Expected Umi: %s. Observed: %s.\n", dna.BasesToString(scAln.Umi), v.ExpectedUmi[i])
			} else if dna.BasesToString(scAln.Bx) != v.ExpectedBx[i] {
				t.Errorf("Error in ToSingleCellAlignment. Expected Bx: %s. Observed: %s.\n", dna.BasesToString(scAln.Bx), v.ExpectedUmi[i])
			}
		}
	}
}
