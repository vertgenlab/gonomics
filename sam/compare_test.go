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
		Equal(Sam{}, r002) ||
		Equal(r004, r002) ||
		Equal(r003, r001Supplemental) ||
		!Equal(Sam{}, Sam{}) {
		t.Error("problem with equal function")
	}

	// nominal test for each field
	if Equal(Sam{QName: "test"}, Sam{}) ||
		Equal(Sam{Flag: 7}, Sam{}) ||
		Equal(Sam{RName: "test"}, Sam{}) ||
		Equal(Sam{Pos: 7}, Sam{}) ||
		Equal(Sam{MapQ: 7}, Sam{}) ||
		Equal(Sam{Cigar: cigar.FromString("7M")}, Sam{}) ||
		Equal(Sam{RNext: "test"}, Sam{}) ||
		Equal(Sam{PNext: 7}, Sam{}) ||
		Equal(Sam{TLen: 7}, Sam{}) ||
		Equal(Sam{Seq: dna.StringToBases("ATG")}, Sam{}) ||
		Equal(Sam{Qual: "FFF"}, Sam{}) ||
		Equal(Sam{Extra: "SA:8"}, Sam{}) {
		t.Error("problem with equal function")
	}
}
