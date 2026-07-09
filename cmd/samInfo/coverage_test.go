package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

/*
Testing files were generated with the following commands in gonomics v1.0.0

~/go/bin/randSeq -lenSeq 1000 -numSeq 1 -setSeed 5 myRandSeq.fa

~/go/bin/simulateSam -coverage 1.0 -flatErrorRate 0.01 -setSeed 9 myRandSeq.fa test1.unsorted.bam
samtools sort test1.unsorted.bam -o test1.bam
rm test1.unsorted.bam

~/go/bin/simulateSam -coverage 10.0 -flatErrorRate 0.01 -setSeed 2 myRandSeq.fa test2.unsorted.bam
samtools sort test2.unsorted.bam -o test2.bam
rm test2.unsorted.bam

~/go/bin/simulateSam -coverage 30.0 -flatErrorRate 0.01 -setSeed 7 myRandSeq.fa test3.unsorted.bam
samtools sort test3.unsorted.bam -o test3.bam
rm test3.unsorted.bam

rm myRandSeq.fa
*/

var SamCoverageTests = []struct {
	InFile           string
	OutHistFile      string
	StatSummaryFile  string
	HighEndFilter    float64
	CountNinDepth    bool
	Verbose          int
	ExpectedHistFile string
	ExpectedStatFile string
}{
	{InFile: "testdata/coverage/test1.bam",
		OutHistFile:      "testdata/coverage/tmp.test1.hist.txt",
		StatSummaryFile:  "testdata/coverage/tmp.test1.stats.txt",
		HighEndFilter:    0.1,
		CountNinDepth:    false,
		Verbose:          0,
		ExpectedHistFile: "testdata/coverage/expected.test1.hist.txt",
		ExpectedStatFile: "testdata/coverage/expected.test1.stats.txt"},
	{InFile: "testdata/coverage/test2.bam",
		OutHistFile:      "testdata/coverage/tmp.test2.hist.txt",
		StatSummaryFile:  "testdata/coverage/tmp.test2.stats.txt",
		HighEndFilter:    0.5,
		CountNinDepth:    false,
		Verbose:          0,
		ExpectedHistFile: "testdata/coverage/expected.test2.hist.txt",
		ExpectedStatFile: "testdata/coverage/expected.test2.stats.txt"},
	{InFile: "testdata/coverage/test3.bam",
		OutHistFile:      "testdata/coverage/tmp.test3.hist.txt",
		StatSummaryFile:  "testdata/coverage/tmp.test3.stats.txt",
		HighEndFilter:    0.01,
		CountNinDepth:    false,
		Verbose:          0,
		ExpectedHistFile: "testdata/coverage/expected.test3.hist.txt",
		ExpectedStatFile: "testdata/coverage/expected.test3.stats.txt"},
}

func TestSamCoverage(t *testing.T) {
	var s CoverageSettings
	for _, v := range SamCoverageTests {
		s = CoverageSettings{
			SamFileName:      v.InFile,
			HistogramOutFile: v.OutHistFile,
			StatSummaryFile:  v.StatSummaryFile,
			HighEndFilter:    v.HighEndFilter,
			CountNinDepth:    v.CountNinDepth,
			Verbose:          v.Verbose,
		}
		samCoverage(s)
		if !fileio.AreEqual(v.OutHistFile, v.ExpectedHistFile) {
			t.Errorf("Error: samInfo coverage. Output histogram was not as expected.\n")
		} else {
			fileio.MustRemove(v.OutHistFile)
		}
		if !fileio.AreEqual(v.StatSummaryFile, v.ExpectedStatFile) {
			t.Errorf("Error: samInfo coverage. Output summary statistics file was not as expected.\n")
		} else {
			fileio.MustRemove(v.StatSummaryFile)
		}
	}
}
