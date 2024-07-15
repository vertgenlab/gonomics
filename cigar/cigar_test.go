package cigar

import (
	"strings"
	"testing"
)

func TestFromString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: Match}
	var c2 Cigar = Cigar{RunLength: 2, Op: Insertion}
	var c3 Cigar = Cigar{RunLength: 16, Op: Deletion}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := FromString(cigarsString)

	for i := 0; i < len(cigarCheck); i++ {
		if !isEqual(cigars[i], cigarCheck[i]) {
			t.Errorf("Error with FromString")
		}
	}
}

func TestToString(t *testing.T) {
	var cigarsString = "35M2I16D"
	var c1 Cigar = Cigar{RunLength: 35, Op: Match}
	var c2 Cigar = Cigar{RunLength: 2, Op: Insertion}
	var c3 Cigar = Cigar{RunLength: 16, Op: Deletion}
	var cigars []Cigar = []Cigar{c1, c2, c3}

	cigarCheck := ToString(cigars)

	if strings.Compare(cigarCheck, cigarsString) != 0 {
		t.Errorf("Error with ToString")
	}
}

func isEqual(a Cigar, b Cigar) bool {
	return (a.RunLength == b.RunLength && a.Op == b.Op)
}

func TestUint32ToCigar(t *testing.T) {
	var query []uint32 = []uint32{560, 33, 258}
	byteCigs := Uint32ToCigar(query)
	if ToString(byteCigs) != "35M2I16D" {
		t.Errorf("Error: 35M!=560, 2I!=33, 16D!=258...\n")
	}
}

func TestToUint32(t *testing.T) {
	var answer []uint32 = []uint32{560, 33, 258}
	byteCigs := []Cigar{c1, c2, c3}
	code := CigarToUint32(byteCigs)
	if len(answer) == len(code) {
		if code[0] != answer[0] || code[1] != answer[1] || code[2] != answer[2] {
			t.Errorf("Error: %d!=35M || %d!=2I || %d!=16D\n", code[0], code[1], code[2])
		}
	}
}

func BenchmarkCigarToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ToString([]Cigar{c1, c2, c3})
	}
}

func BenchmarkStringToCigar(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var cigarstring string = "35M2I16D"
	for n := 0; n < b.N; n++ {
		FromString(cigarstring)
	}
}
