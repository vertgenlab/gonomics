package sam

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"sort"
	"strings"
)

// Equal returns true if the two input Sam structs are identical.
func Equal(a Sam, b Sam) bool {
	if a.QName != b.QName {
		return false
	}
	if a.Flag != b.Flag {
		return false
	}
	if a.RName != b.RName {
		return false
	}
	if a.Pos != b.Pos {
		return false
	}
	if a.MapQ != b.MapQ {
		return false
	}
	if cigar.ToString(a.Cigar) != cigar.ToString(b.Cigar) {
		return false
	}
	if a.RNext != b.RNext {
		return false
	}
	if a.PNext != b.PNext {
		return false
	}
	if a.TLen != b.TLen {
		return false
	}
	if dna.CompareSeqsIgnoreCase(a.Seq, b.Seq) != 0 {
		return false
	}
	if a.Qual != b.Qual {
		return false
	}
	if a.Extra != b.Extra {
		return false
	}
	return true
}

func SortByCoord(a []Sam) {
	sort.Slice(a, func(i, j int) bool { return Compare(a[i], a[j]) == -1 })
}

func Compare(a Sam, b Sam) int {
	chromComp := strings.Compare(a.RName, b.RName)
	if chromComp != 0 {
		return chromComp
	}
	if a.Pos < b.Pos {
		return -1
	}
	if a.Pos > b.Pos {
		return 1
	} //do I need to add anything about alignment lengths to break start position ties?
	return 0

}