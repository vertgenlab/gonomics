package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
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
	InFile        string
	OutFile       string
	CountNinDepth bool
	Verbose       int
	ExpectedFile  string
}{
	{InFile: "testdata/test1.bam",
		OutFile:       "testdata/tmp.test1.txt",
		CountNinDepth: false,
		Verbose:       0,
		ExpectedFile:  "testdata/expected.test1.txt"},
	{InFile: "testdata/test2.bam",
		OutFile:       "testdata/tmp.test2.txt",
		CountNinDepth: false,
		Verbose:       0,
		ExpectedFile:  "testdata/expected.test2.txt"},
	{InFile: "testdata/test3.bam",
		OutFile:       "testdata/tmp.test3.txt",
		CountNinDepth: false,
		Verbose:       0,
		ExpectedFile:  "testdata/expected.test3.txt"},
}

func TestSamCoverage(t *testing.T) {
	var s Settings
	for _, v := range SamCoverageTests {
		s = Settings{
			SamFileName:   v.InFile,
			OutFile:       v.OutFile,
			CountNinDepth: v.CountNinDepth,
			Verbose:       v.Verbose,
		}
		samCoverage(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in samCoverage")
		} else {
			err := os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
