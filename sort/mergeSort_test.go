package sort

import (
	"github.com/vertgenlab/gonomics/vcf"
	"math/rand"
	"sort"
	"testing"
)

func TestMergeSort(t *testing.T) {
	rand.Seed(1)
	vcfs, _ := vcf.Read("testdata/test.vcf")
	data := make(chan vcf.Vcf, len(vcfs))
	rand.Shuffle(len(vcfs), func(i, j int) {
		vcfs[i], vcfs[j] = vcfs[j], vcfs[i]
	})

	for i := range vcfs {
		data <- vcfs[i]
	}
	close(data)

	out := GoExternalMergeSort[vcf.Vcf](data, 5, func(a, b vcf.Vcf) bool {
		return a.Pos < b.Pos
	})

	answer := make([]vcf.Vcf, 0, len(vcfs))
	for v := range out {
		answer = append(answer, v)
	}

	if len(answer) != len(vcfs) {
		t.Error("problem with mergeSort")
	}

	sort.Slice(vcfs, func(i, j int) bool {
		return vcfs[i].Pos < vcfs[j].Pos
	})
	if !vcf.AllEqual(vcfs, answer) {
		t.Error("problem with mergeSort")
	}
}
