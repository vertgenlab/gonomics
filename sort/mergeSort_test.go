package sort

import "testing"

func TestExternalMergeSort(t *testing.T) {
	ExternalMergeSort("testdata/test.vcf", 100, "tmp", "sorted.vcf")
}
