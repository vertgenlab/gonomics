package interval

import (
	"testing"

	"github.com/vertgenlab/gonomics/bed"
)

var CoordsToStringTests = []struct {
	I        Interval
	Expected string
}{
	{I: bed.Bed{Chrom: "chr1", ChromStart: 2900, ChromEnd: 3100}, Expected: "chr1:2900-3100"},
}

func TestCoordsToString(t *testing.T) {
	var currAnswer string
	for _, v := range CoordsToStringTests {
		currAnswer = CoordsToString(v.I)
		if currAnswer != v.Expected {
			t.Errorf("Error in CoordsToString. Output is not as expected.")
		}
	}
}
