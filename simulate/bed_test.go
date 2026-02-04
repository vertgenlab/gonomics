package simulate

import (
	"math/rand"
	"os"
	"testing"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var CountWindowsTests = []struct {
	InputFa      string
	RegionLength int
	Expected     int
}{
	{"testdata/ref_short.fasta", 3, 59},        // regionlength>1
	{"testdata/ref_short.fasta", 1, 69},        // regionlength=1
	{"testdata/ref_short_allGaps.fasta", 5, 0}, // all gaps
	{"testdata/ref_short_2.fasta", 50, 1},      // all gaps
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
	InputFa      string
	Pos          int
	RegionLength int
	Expected     string
}{
	//{"testdata/ref_short.fasta", 49, 1, "testdata/generateBedRegion_expected_1.bed"},   // end of region, len=1
	//{"testdata/ref_short.fasta", 10, 3, "testdata/generateBedRegion_expected_2.bed"},   // end of region before gap, len>1
	{"testdata/ref_short.fasta", 0, 1, "testdata/generateBedRegion_expected_3.bed"},    // beginning of region, len=1
	{"testdata/ref_short.fasta", 50, 1, "testdata/generateBedRegion_expected_4.bed"},   // beginning of new region, len=1
	{"testdata/ref_short_2.fasta", 0, 50, "testdata/generateBedRegion_expected_5.bed"}, // entire region
	{"testdata/ref_short.fasta", 14, 13, "testdata/generateBedRegion_expected_6.bed"},  // after gap, to end of region, len>1
}

func TestGenerateBedRegion(t *testing.T) {
	var searchSpace []bed.Bed
	var inputFile []fasta.Fasta
	var region bed.Bed
	var generatedBool bool
	var expected bed.Bed
	for idx, v := range GenerateBedRegionTests {
		inputFile = fasta.Read(v.InputFa)
		searchSpace = bed.UngappedRegionsAllFromFa(inputFile)
		region, generatedBool = GenerateBedRegion(searchSpace, v.Pos, v.RegionLength)
		expected = bed.Read(v.Expected)[0]
		if !generatedBool {
			t.Errorf("No bed generated.")
		}
		if !bed.Equal(region, expected) {
			fmt.Println(region)
			t.Errorf("Error in GenerateBedRegion testcase %v.", idx)
		}
	}
}

var GoSimulateBedTests = []struct {
	InputFa      string
	RegionCount  int
	RegionLength int
	SetSeed      int64
	Expected     string
	OutFile      string
}{
	{"testdata/ref_short.fasta", 3, 1, 8, "testdata/goSimulateBed_expected_1.bed", "testdata/goSimulateBed_out_1.bed"},
	{"testdata/ref_short.fasta", 100, 3, 3, "testdata/goSimulateBed_expected_2.bed", "testdata/goSimulateBed_out_2.bed"},
}

func TestGoSimulateBed(t *testing.T) {
	var searchSpace []bed.Bed
	var inputFile []fasta.Fasta
	var err error
	for idx, v := range GoSimulateBedTests {
		rand.Seed(v.SetSeed)
		inputFile = fasta.Read(v.InputFa)
		searchSpace = bed.UngappedRegionsAllFromFa(inputFile)
		c := GoSimulateBed(searchSpace, v.RegionCount, v.RegionLength)
		out := fileio.EasyCreate(v.OutFile)
		for i := range c {
			bed.WriteBed(out, i)
		}
		err = out.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in GoSimulateBed testcase %v", idx)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
