package bed

import (
	"fmt"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
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

var inB []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Score: 10, Name: "A"}, {Chrom: "chr1", ChromStart: 15, ChromEnd: 25, Score: 40, Name: "B"}, {Chrom: "chr1", ChromStart: 25, ChromEnd: 45, Score: 500, Name: "C"}}
var expectedMergeFalse []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 25, Score: 40}, {Chrom: "chr1", ChromStart: 25, ChromEnd: 45, Score: 500}}
var expectedMergeTrue []Bed = []Bed{{Chrom: "chr1", ChromStart: 10, ChromEnd: 45, Score: 500}}

var MergeHighMemTests = []struct {
	InBed         []Bed
	ExpectedBed   []Bed
	MergeAdjacent bool
	KeepAllNames  bool
}{
	{InBed: inB, ExpectedBed: expectedMergeFalse, MergeAdjacent: false, KeepAllNames: false},
	{InBed: inB, ExpectedBed: expectedMergeTrue, MergeAdjacent: true, KeepAllNames: false},
}

func TestMergeHighMem(t *testing.T) {
	var outB []Bed
	for _, v := range MergeHighMemTests {
		outB = MergeHighMem(v.InBed, v.MergeAdjacent, v.KeepAllNames)
		if !AllAreEqual(outB, v.ExpectedBed) { //doesn't care about name
			fmt.Printf("MergeAdjacent: %v.\n", v.MergeAdjacent)
			for i := range outB {
				fmt.Printf("%s\n", ToString(outB[i], 5))
			}

			t.Errorf("Error in MergeHighMem. Output was not as expected.")
		}
	}
}

var MergeKeepNamesTest = []struct {
	InFile        string
	ExpectedFile  string
	MergeAdjacent bool
	KeepAllNames  bool
}{
	{"testdata/test.names.bed", "testdata/test.names.merged.bed", false, true},
	{"testdata/test.names.bed", "testdata/test.names.adjacent.merged.bed", true, true},
}

func TestMergeKeepNames(t *testing.T) {
	var outB []Bed
	var err error
	for _, i := range MergeKeepNamesTest {
		outB = MergeHighMem(Read(i.InFile), i.MergeAdjacent, i.KeepAllNames)
		Write("testdata/tmp.txt", outB)
		if !fileio.AreEqual("testdata/tmp.txt", i.ExpectedFile) {
			t.Errorf("Error in bedMerge.")
		} else {
			err = os.Remove("testdata/tmp.txt")
			exception.PanicOnErr(err)
		}
	}
}
