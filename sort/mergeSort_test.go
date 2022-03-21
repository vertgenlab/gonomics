package sort

import (
	"fmt"
	"github.com/vertgenlab/gonomics/vcf"
	"math/rand"
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

	for v := range out {
		fmt.Println(v.Chr, v.Pos, v.Alt)
	}
}
