package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"testing"
)

var genes = []bed.Bed{{Chrom: "chr1", ChromStart: 2, ChromEnd: 3, Name: "first", Score: 0}, {Chrom: "chr1", ChromStart: 13, ChromEnd: 14, Name: "second", Score: 0}, {Chrom: "chr1", ChromStart: 50, ChromEnd: 51, Name: "third", Score: 0}, {Chrom: "chr2", ChromStart: 10, ChromEnd: 40, Name: "fourth", Score: 0}}
var bedpe = []BedPe{{bed.Bed{Chrom: "chr1", ChromStart: 0, ChromEnd: 2, Name: "", Score: 0}, bed.Bed{Chrom: "chr1", ChromStart: 8, ChromEnd: 11, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr2", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr2", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}, {bed.Bed{Chrom: "chr3", ChromStart: 0, ChromEnd: 5, Name: "", Score: 0}, bed.Bed{Chrom: "chr3", ChromStart: 85, ChromEnd: 95, Name: "", Score: 0}}}
var size = []chromInfo.ChromInfo{{Name: "chr1", Size: 600}, {Name: "chr2", Size: 100}}

func TestFill3dSpace(t *testing.T) {
	answer := Fill3dSpace(bedpe, genes, chromInfo.SliceToMap(size))
	bed.Write("testdata/fill3dSpaceOut.bed", answer)

}
