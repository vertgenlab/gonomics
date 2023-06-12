package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
)

var MultiFaExtractTests = []struct {
	inputFile    string
	outputFile   string
	expectedFile string
	removeGaps   bool
	bed          string
	start        int
	end          int
	out1         string //these six fields are for the three output and expected files made by the bed option usage of this cmd.
	out2         string
	out3         string
	expected1    string
	expected2    string
	expected3    string
}{
	{"testdata/testInput.fa", "testdata/testOut.fa", "testdata/testOut.10to200.fa", false, "", 10, 200, "", "", "", "", "", ""},
	{"testdata/testInput.fa", "testdata/testOut.fa", "testdata/testOut.10to200.RemoveGaps.fa", true, "", 10, 200, "", "", "", "", "", ""},
	{"testdata/testInput.fa", "", "", false, "testdata/test.bed", -1, -1, "chr1.20.30.fa", "chr1.30.50.fa", "chr1.60.200.fa", "testdata/chr1.20.30.expected.fa", "testdata/chr1.30.50.expected.fa", "testdata/chr1.60.200.expected.fa"},                  //~/go/bin/multiFaExtract -bed testdata/test.bed testdata/testInput.fa
	{"testdata/testInput.fa", "", "", true, "testdata/test.bed", -1, -1, "chr1.20.30.fa", "chr1.30.50.fa", "chr1.60.200.fa", "testdata/chr1.20.30.expected.noGap.fa", "testdata/chr1.30.50.expected.noGap.fa", "testdata/chr1.60.200.expected.noGap.fa"}, //~/go/bin/multiFaExtract -bed testdata/test.bed -removeGaps testdata/testInput.fa
}

func TestMultiFaExtract(t *testing.T) {
	var err error
	var s Settings
	for _, v := range MultiFaExtractTests {
		s = Settings{
			InFile:     v.inputFile,
			OutFile:    v.outputFile,
			Start:      v.start,
			End:        v.end,
			Bed:        v.bed,
			RemoveGaps: v.removeGaps,
		}
		multiFaExtract(s)
		if v.bed == "" {
			records := fasta.Read(v.outputFile)
			expected := fasta.Read(v.expectedFile)
			if !fasta.AllAreEqual(records, expected) {
				t.Errorf("Error in multiFaExtract.")
			} else {
				err = os.Remove(v.outputFile)
				exception.PanicOnErr(err)
			}
		} else {
			rec1 := fasta.Read(v.out1)
			exp1 := fasta.Read(v.expected1)
			if !fasta.AllAreEqual(rec1, exp1) {
				t.Errorf("Error in multiFaExtract, bed usage. First outFile.")
			} else {
				err = os.Remove(v.out1)
				exception.PanicOnErr(err)
			}
			rec2 := fasta.Read(v.out2)
			exp2 := fasta.Read(v.expected2)
			if !fasta.AllAreEqual(rec2, exp2) {
				t.Errorf("Error in multiFaExtract, bed usage. Second outFile.")
			} else {
				err = os.Remove(v.out2)
				exception.PanicOnErr(err)
			}
			rec3 := fasta.Read(v.out3)
			exp3 := fasta.Read(v.expected3)
			if !fasta.AllAreEqual(rec3, exp3) {
				t.Errorf("Error in multiFaExtract, bed usage. Third outFile.")
			} else {
				err = os.Remove(v.out3)
				exception.PanicOnErr(err)
			}
		}
	}
}
