package interval

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"os"
	"testing"
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

func TestIntervalSimilarity(t *testing.T) {
	var expA float64 = 0.75
	var expB float64 = float64(2) / float64(3)
	var expAvg float64 = numbers.AverageFloat64([]float64{expA, expB})

	a := bed.Read("../cmd/bedSimilarity/testdata/smallAJ.bed")
	b := bed.Read("../cmd/bedSimilarity/testdata/largeAJ.bed")
	aSlice := BedSliceToIntervals(a)
	bSlice := BedSliceToIntervals(b)
	ansA, ansB, ansAvg := IntervalSimilarity(aSlice, bSlice)
	if expA != ansA {
		t.Errorf("expected proportionA does not equal actual proportionA. Expected: %f. Actual: %f", expA, ansA)
	}
	if expB != ansB {
		t.Errorf("expected proportionB does not equal actual proportionB. Expected: %f. Actual: %f", expB, ansB)
	}
	if expAvg != ansAvg {
		t.Errorf("expected averageProportion does not equal actual averageProportion. Expected: %f. Actual: %f", expAvg, ansAvg)
	}
}

func TestAreEqual(t *testing.T) {
	in := bed.Read("testdata/in.bed")
	if !AreEqual(in[0], in[3]) {
		t.Error("Expected to return equal for these intevals,", in[0], in[3])
	}
	if AreEqual(in[0], in[2]) {
		t.Error("Expected to return un-equal for these intevals,", in[0], in[2])
	}
}

func TestUnique(t *testing.T) {
	in := bed.Read("testdata/in.bed")
	out := fileio.EasyCreate("testdata/out.unique.bed")
	inIntervals := BedSliceToIntervals(in)
	outIntervals := Unique(inIntervals)
	for _, i := range outIntervals {
		fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%d\t%d", i.GetChrom(), i.GetChromStart(), i.GetChromEnd()))
	}
	err := out.Close()
	exception.PanicOnErr(err)
	if !fileio.AreEqual("testdata/exp.unique.bed", "testdata/out.unique.bed") {
		t.Errorf("Error in unique.go. testdata/exp.unique.bed and testdata/out.unique.bed aren't equal to each other")
	} else {
		err = os.Remove("testdata/out.unique.bed")
		exception.PanicOnErr(err)
	}
}
