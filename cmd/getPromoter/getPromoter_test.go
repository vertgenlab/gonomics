package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"os"
	"testing"
)

func TestGetPromoter(t *testing.T) {

	getPromoter("testdata/uniqueGenes.txt", "testdata/gtfFileTest.gtf", "testdata/tmp.bed", "testdata/hg38.chrom.sizes", 1000, 200)
	if !bed.AllAreEqual(bed.Read("testdata/expected1kb.bed"), bed.Read("testdata/tmp.bed")) {
		t.Errorf("bed output was not as expected. See tmp.bed")
	}
	getPromoter("testdata/uniqueGenes.txt", "testdata/gtfFileTest.gtf", "testdata/tmp.bed", "testdata/hg38.chrom.sizes", 5000, 1000)
	if !bed.AllAreEqual(bed.Read("testdata/expected5kb.bed"), bed.Read("testdata/tmp.bed")) {
		t.Errorf("bed output was not as expected. See tmp.bed")
	} else {
		os.Remove("testdata/tmp.bed")
	}
}
