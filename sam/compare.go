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

func SortByCoord(a []Sam, head Header) {
	sort.Slice(a, func(i, j int) bool { return Compare(a[i], a[j], head) == -1 })
}

func Compare(a Sam, b Sam, head Header) int {
	var chromComp int
	if head.Chroms != nil && a.RName != b.RName {
		chromComp = chromCompFromHeader(a.RName, b.RName, head)
	} else {
		chromComp = strings.Compare(a.RName, b.RName)
	}
	if chromComp != 0 {
		return chromComp
	}
	if a.Pos < b.Pos {
		return -1
	}
	if a.Pos > b.Pos {
		return 1
	} //do I need to add anything about alignment lengths to break start position ties? Look at how samtools sorts exactly.
	return 0
}

func chromCompFromHeader(aName string, bName string, head Header) int {
	var aOrder, bOrder int
	var aFound, bFound bool = false, false
	for _, i := range head.Chroms {
		if i.Name == aName {
			aOrder = i.Order
			aFound = true
		}
		if i.Name == bName {
			bOrder = i.Order
			bFound = true
		}
		if aFound && bFound {
			break
		}
	}
	if aOrder < bOrder {
		return -1
	}
	if aOrder > bOrder {
		return 1
	}
	return 0
}
