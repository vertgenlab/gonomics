package simulate

import (
	// "os"
	"testing"
	"math/rand"

	// "github.com/vertgenlab/gonomics/exception"
	// "github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"

)

var CountWindowsTests = []struct {
	InputFa string
	RegionLength int
	Expected int
}{
	{"testdata/ref_short.fasta", 3, 59}, // regionlength>1
	{"testdata/ref_short.fasta", 1, 69}, // regionlength=1
	{"testdata/ref_short_allGaps.fasta", 5, 0}, // all gaps
	{"testdata/ref_short_2.fasta", 50, 1}, // all gaps
}

func TestCountWindows(t *testing.T) {
	var searchSpace []bed.Bed
	var numWindows int
	var inputFile []fasta.Fasta
	for _, v := range CountWindowsTests {
		inputFile = fasta.Read(v.InputFa)
		searchSpace = bed.UngappedRegionsAllFromFa(inputFile)
		numWindows = CountWindows(searchSpace, v.RegionLength)
		if numWindows != v.Expected {
			t.Errorf("Error in CountWindows. Output %v did not match expected %v.", v.Expected, numWindows)
		}
	}
}

var GenerateBedRegionTests = []struct {
	InputFa string
	Pos int
	RegionLength int
	SetSeed int64
	Expected string
} {
	{"testdata/ref_short.fasta", 51, 1, 8, "testdata/generateBedRegion_expected_1.bed"}, // end of region, len=1
	{"testdata/ref_short.fasta", 49, 3, 3, "testdata/generateBedRegion_expected_2.bed"}, // end of region, len>1
	{"testdata/ref_short.fasta", 0, 1, 22, "testdata/generateBedRegion_expected_3.bed"}, // beginning of region, len=1
	{"testdata/ref_short.fasta", 52, 1, 5, "testdata/generateBedRegion_expected_4.bed"}, // beginning of new region, len=1
	{"testdata/ref_short_2.fasta", 0, 50, 2, "testdata/generateBedRegion_expected_5.bed"}, // entire region
	{"testdata/ref_short.fasta", 15, 13, 17, "testdata/generateBedRegion_expected_5.bed"}, // after gap, to end of region, len>1
}

func TestGenerateBedRegion(t * testing.T) {
	var searchSpace []bed.Bed
	var inputFile []fasta.Fasta
	// var region bed.Bed
	for idx, v := range GenerateBedRegionTests {
		rand.Seed(v.SetSeed)
		inputFile = fasta.Read(v.InputFa)
		searchSpace = bed.UngappedRegionsAllFromFa(inputFile)
		log.Printf("test number %v, input posiion: %v, regionlength %v\n", idx, v.Pos, v.RegionLength)
		GenerateBedRegion(searchSpace, v.Pos, v.RegionLength)
		// region = GenerateBedRegion(searchSpace, v.Pos, v.RegionLength)
		// if !bed.AllAreEqual(region, v.Expected) {
		// 	t.Errorf("Error in generateBedRegion.")
		// }
	}
}
