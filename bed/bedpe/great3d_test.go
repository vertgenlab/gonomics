package bedpe

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
)

var genes = []bed.Bed{{Chrom: "chr1", ChromStart: 2, ChromEnd: 3, Name: "first", Score: 0}, {Chrom: "chr1", ChromStart: 13, ChromEnd: 14, Name: "second", Score: 0}, {Chrom: "chr1", ChromStart: 500, ChromEnd: 501, Name: "third", Score: 0}, {Chrom: "chr2", ChromStart: 10, ChromEnd: 40, Name: "fourth", Score: 0}}
var bedpe = []BedPe{{bed.Bed{Chrom: "chr1", ChromStart: 80, ChromEnd: 81, Name: "", Score: 0}, bed.Bed{Chrom: "chr1", ChromStart: 300, ChromEnd: 301, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr2", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr2", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr3", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr3", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}}
var size = []chromInfo.ChromInfo{{Name: "chr1", Size: 600}, {Name: "chr2", Size: 100}}

func TestFill3dSpace(t *testing.T) {
	answer := Fill3dSpace(bedpe, genes, chromInfo.SliceToMap(size))
	bed.Write("testdata/fill3dSpaceOut.bed", answer)

	if !bed.AllAreEqual(answer, bed.Read("testdata/expected.out")) {
		t.Errorf("Error: output didn't match expected file.")
	} else {
		os.Remove("testdata/fill3dSpaceOut.bed")
	}
}
