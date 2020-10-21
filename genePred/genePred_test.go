package genePred

import (
	"testing"
)

var genePred1 = GenePred{Id: "test", Chrom: "0", TxStart: 0, TxEnd: 1000, CdsStart: 0, CdsEnd: 900, ExonStarts: []int{0, 18, 500, 800}, ExonEnds: []int{2, 57, 601, 830}}
var genePred2 = GenePred{Id: "test", Chrom: "1", TxStart: 0, TxEnd: 1000, CdsStart: 0, CdsEnd: 900, ExonStarts: []int{0, 18, 500, 800}, ExonEnds: []int{2, 57, 601, 830}}
var genePreds []*GenePred = []*GenePred{&genePred1, &genePred2}

var ReadTests = []struct {
	name string
	data []*GenePred
}{
	{"testGenePred.gpd", genePreds},
}

func TestRead(t *testing.T) {
	for _, test := range ReadTests {
		actual := Read(test.name)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.name)
		}
	}
}
