package simulate

import (
	"math/rand"
	"os"
	"testing"

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
	{
		InputFa:      "testdata/ref_short.fasta",
		RegionLength: 3,
		Expected:     59,
	}, // region length > 1
	{
		InputFa:      "testdata/ref_short.fasta",
		RegionLength: 1,
		Expected:     69,
	}, // region length = 1
	{
		InputFa:      "testdata/ref_short_allGaps.fasta",
		RegionLength: 5,
		Expected:     0,
	}, // all gaps
	{
		InputFa:      "testdata/ref_short_2.fasta",
		RegionLength: 50,
		Expected:     1,
	}, // single region
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
	{
		InputFa:      "testdata/ref_short.fasta",
		Pos:          49,
		RegionLength: 1,
		Expected:     "testdata/generateBedRegion_expected_1.bed",
	}, // End of region, length = 1
	{
		InputFa:      "testdata/ref_short.fasta",
		Pos:          10,
		RegionLength: 3,
		Expected:     "testdata/generateBedRegion_expected_2.bed",
	}, // End of region before gap, length > 1
	{
		InputFa:      "testdata/ref_short.fasta",
		Pos:          0,
		RegionLength: 1,
		Expected:     "testdata/generateBedRegion_expected_3.bed",
	}, // Beginning of region, length = 1
	{
		InputFa:      "testdata/ref_short.fasta",
		Pos:          50,
		RegionLength: 1,
		Expected:     "testdata/generateBedRegion_expected_4.bed",
	}, // Beginning of new region, length = 1
	{
		InputFa:      "testdata/ref_short_2.fasta",
		Pos:          0,
		RegionLength: 50,
		Expected:     "testdata/generateBedRegion_expected_5.bed",
	}, // Entire region
	{
		InputFa:      "testdata/ref_short.fasta",
		Pos:          14,
		RegionLength: 13,
		Expected:     "testdata/generateBedRegion_expected_6.bed",
	}, // After gap, to end of region, length > 1
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
	{
		InputFa:      "testdata/ref_short.fasta",
		RegionCount:  3,
		RegionLength: 1,
		SetSeed:      8,
		Expected:     "testdata/goSimulateBed_expected_1.bed",
		OutFile:      "testdata/goSimulateBed_out_1.bed",
	},
	{
		InputFa:      "testdata/ref_short.fasta",
		RegionCount:  100,
		RegionLength: 3,
		SetSeed:      3,
		Expected:     "testdata/goSimulateBed_expected_2.bed",
		OutFile:      "testdata/goSimulateBed_out_2.bed",
	},
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
