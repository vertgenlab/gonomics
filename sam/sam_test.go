package sam

import (
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
)

var s Sam = Sam{
	QName: "read",
	Flag:  147,
	MapQ:  30,
	RName: "ref",
	Pos:   37,
	Cigar: cigar.FromString("9M"),
	RNext: "=",
	PNext: 7,
	TLen:  -39,
	Seq:   dna.StringToBases("CAGCGGCAT"),
	Qual:  "*",
	Extra: "NM:i:1",
}

var expected string = "read\t147\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t*\tNM:i:1"

func TestSamString(t *testing.T) {
	result := s.String()
	if result != expected {
		t.Errorf("Error: String() = %q, want %q", result, expected)
	}

}

func BenchmarkSamString(b *testing.B) {
	b.ReportAllocs()

	for i := 0; i < b.N; i++ {
		result := s.String()
		if result != expected {
			b.Errorf("Error: String() = %q, want %q", result, expected)
		}
	}
}
