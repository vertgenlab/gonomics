package starrSeq

import (
	"testing"
)

func TestSortReadByCellBx(t *testing.T) {
	var pre []Read = []Read{{Bx: "GATC", Construct: "1"}, {Bx: "CATG", Construct: "2"}, {Bx: "GATC", Construct: "3"}, {Bx: "AAAA", Construct: "4"}}
	var exp []Read = []Read{{Bx: "AAAA", Construct: "4"}, {Bx: "CATG", Construct: "2"}, {Bx: "GATC", Construct: "1"}, {Bx: "GATC", Construct: "3"}}

	SortReadByCellBx(pre)

	for i := range pre {
		if pre[i] != exp[i] {
			t.Errorf("ERROR: sorting Read slices failed. %s %s", pre[i], exp[i])
		}
	}
}
