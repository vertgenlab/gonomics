package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var TrimTests = []struct {
	InFile       string
	ExpectedFile string
	TrimLeft     int
	TrimRight    int
}{
	{"testdata/testTrim.bed", "testdata/expectedTrim.bed", 10, 10},
}

func TestTrim(t *testing.T) {
	var b, expected []Bed
	for _, v := range TrimTests {
		b = Read(v.InFile)
		expected = Read(v.ExpectedFile)
		Trim(b, v.TrimLeft, v.TrimRight)
		if !AllAreEqual(b, expected) {
			t.Errorf("Error in Trim.")
		}
	}
}

var inB []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Score: 10}, {Chrom: "chr1", ChromStart: 15, ChromEnd: 25, Score: 40}, {Chrom: "chr1", ChromStart: 25, ChromEnd: 45, Score: 500}}
var expectedBfalse []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 25, Score: 40}, {Chrom: "chr1", ChromStart: 25, ChromEnd: 45, Score: 500}}
var expectedBtrue []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 45, Score: 500}}

var MergeHighMemTests = []struct {
	InBed         []Bed
	ExpectedBed   []Bed
	MergeAdjacent bool
}{
	{InBed: inB, ExpectedBed: expectedBfalse, MergeAdjacent: false},
	{InBed: inB, ExpectedBed: expectedBtrue, MergeAdjacent: true},
}

func TestMergeHighMem(t *testing.T) {
	var outB []Bed
	for _, v := range MergeHighMemTests {
		outB = MergeHighMem(v.InBed, v.MergeAdjacent)
		if !AllAreEqual(outB, v.ExpectedBed) {
			fmt.Printf("MergeAdjacent: %v.\n", v.MergeAdjacent)
			for i := range outB {
				fmt.Printf("%s\n", ToString(outB[i], 5))
			}

			t.Errorf("Error in MergeHighMem. Output was not as expected.")
		}
	}
}

var FillSpaceTests = []struct {
	InputFile string
	Genome    map[string]chromInfo.ChromInfo
	OutFile   string
	Expected  string
}{
	{
		InputFile: "testdata/FillSpace.Input.bed",
		Genome:    map[string]chromInfo.ChromInfo{"chr1": {Name: "chr1", Size: 600}, "chr2": {Name: "chr2", Size: 60}},
		OutFile:   "testdata/tmp.FillSpace.bed",
		Expected:  "testdata/FillSpace.Expected.bed",
	},
}

func TestFillSpaceNoHiddenValue(t *testing.T) {
	var err error
	var records, answer []Bed
	for _, v := range FillSpaceTests {
		records = Read(v.InputFile)
		answer = FillSpaceNoHiddenValue(records, v.Genome)
		Write(v.OutFile, answer)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in FillSpace. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var FillThreeDSpaceTests = []struct {
	InputFile string
	Genome    map[string]chromInfo.ChromInfo
	OutFile   string
	Expected  string
}{
	{
		InputFile: "testdata/FillSpace.Hidden.Input.bed",
		Genome:    map[string]chromInfo.ChromInfo{"chr1": {Name: "chr1", Size: 600}, "chr2": {Name: "chr2", Size: 60}},
		OutFile:   "testdata/tmp.Hidden.FillSpace.bed",
		Expected:  "testdata/FillSpace.Hidden.Expected.bed",
	},
}

func TestFillSpaceHiddenValue(t *testing.T) {
	var err error
	var records, answer []Bed
	for _, v := range FillThreeDSpaceTests {
		records = Read(v.InputFile)
		answer = FillSpaceHiddenValue(records, v.Genome)
		Write(v.OutFile, answer)
		if !fileio.AreEqual(v.OutFile, v.Expected) {
			t.Errorf("Error in FillSpace. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}


