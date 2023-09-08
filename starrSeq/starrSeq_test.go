package starrSeq

import (
	"testing"
)

func TestSortUmiByCellBx(t *testing.T) {
	var pre []UMI = []UMI{{Bx: "GATC", Construct: "1"}, {Bx: "CATG", Construct: "2"}, {Bx: "GATC", Construct: "3"}, {Bx: "AAAA", Construct: "4"}}
	var exp []UMI = []UMI{{Bx: "AAAA", Construct: "4"}, {Bx: "CATG", Construct: "2"}, {Bx: "GATC", Construct: "1"}, {Bx: "GATC", Construct: "3"}}

	SortUmiByCellBx(pre)

	for i := range pre {
		if pre[i] != exp[i] {
			t.Errorf("ERROR: sorting UMI slices failed. %s %s", pre[i], exp[i])
		}
	}
}
