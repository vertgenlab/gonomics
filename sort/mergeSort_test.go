package sort

import (
	"github.com/vertgenlab/gonomics/vcf"
	"testing"
)

//TODO: add test function for each data type
func TestExternalMergeSort(t *testing.T) {
	ExternalMergeSort("testdata/test.vcf", 100, "tmp", "sorted.vcf")
	sorted := vcf.Read("sorted.vcf")
	test := vcf.ByGenomicCoordinates(sorted)
	for i := 1; i < len(sorted); i++ {
		if !test.Less(i-1, i) {
			t.Errorf("ERROR: Problem with external merge sort of vcf files: \n %v \n is not less than \n %v", sorted[i-1], sorted[i])
		}
	}
}
