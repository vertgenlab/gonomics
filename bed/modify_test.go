package bed

import (
	"fmt"
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

var inB []Bed = []Bed{Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 20, Score: 10}, Bed{Chrom: "chr1", ChromStart: 15, ChromEnd: 25, Score: 40}, Bed{Chrom: "chr1", ChromStart:25, ChromEnd:45, Score: 500}}
var expectedBfalse []Bed = []Bed{Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 25, Score: 40}, Bed{Chrom: "chr1", ChromStart:25, ChromEnd:45, Score: 500}}
var expectedBtrue []Bed = []Bed{Bed{Chrom: "chr1", ChromStart: 10, ChromEnd: 45, Score: 500}}

var MergeHighMemTests = []struct {
	InBed []Bed
	ExpectedBed []Bed
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