package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestEqual(t *testing.T) {
	if !Equal(r001, r001) || // find tested records in io_test.go
		Equal(r001, r002) ||
		Equal(r001, r001Supplemental) ||
		Equal(r003, r002) ||
		!Equal(r002, r002) ||
		Equal(r003Supplemental, r003) ||
		Equal(Aln{}, r002) ||
		Equal(r004, r002) ||
		Equal(r003, r001Supplemental) ||
		!Equal(Aln{}, Aln{}) {
		t.Error("problem with equal function")
	}

	if Equal(Aln{QName: "test"}, Aln{}) ||
		Equal(Aln{Flag: 7}, Aln{}) ||
		Equal(Aln{RName: "test"}, Aln{}) ||
		Equal(Aln{Pos: 7}, Aln{}) ||
		Equal(Aln{MapQ: 7}, Aln{}) ||
		Equal(Aln{Cigar: cigar.FromString("7M")}, Aln{}) ||
		Equal(Aln{RNext: "test"}, Aln{}) ||
		Equal(Aln{PNext: 7}, Aln{}) ||
		Equal(Aln{TLen: 7}, Aln{}) ||
		Equal(Aln{Seq: dna.StringToBases("ATG")}, Aln{}) ||
		Equal(Aln{Qual: "FFF"}, Aln{}) ||
		Equal(Aln{Extra: "SA:8"}, Aln{}) {
		t.Error("problem with equal function")
	}
}
