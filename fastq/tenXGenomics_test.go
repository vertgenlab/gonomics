package fastq

import (
	"fmt"
	"testing"
)
//TODO: write a better test
func TestCheck10xBarcodes(t *testing.T) {
	TenxFq := Read("testdata/10x.barcoded_test.fastq")
	for _, read := range TenxFq {
		fmt.Printf("%s\n", fastqStats(read))
	}

}

func TestTrimmingBarcodes(t *testing.T) {
	tenX := make(chan *LinkedRead)
	go ReadToChanLinked("testdata/barcode10x_R1.fastq", "testdata/barcode10x_R2.fastq", tenX)

	for fq := range tenX {
		fmt.Printf("%s\n", PrettyPrint(fq))
	}
}
