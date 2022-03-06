package lift

import (
	"github.com/vertgenlab/gonomics/bed"
	"testing"
)

var ReadingTests = []struct {
	Infile string
	Expected []Lift
}{
	{"testdata/elements1.bed", []Lift{bed.Bed{Chrom: "chr1", ChromStart: 2, ChromEnd: 10}, bed.Bed{Chrom: "chr1", ChromStart: 13, ChromEnd: 15}, bed.Bed{Chrom: "chr1", ChromStart: 50, ChromEnd: 500}}},
}

func TestReading(t *testing.T) {
	var l []Lift
	for _, v := range ReadingTests {
		l = GoRead(v.Infile)
		if !allAreEqual(l, v.Expected) {
			t.Errorf("Error in lift GoRead. Input did not match expected.")
		}
	}
}
