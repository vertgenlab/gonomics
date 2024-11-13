package bed

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ToLowerTests = []struct {
	InFastaFile        string
	InBedFile          string
	IgnoreExtraRegions bool
	OutFile            string
	ExpectedFile       string
}{
	{InFastaFile: "testdata/toLower.fa",
		InBedFile:          "testdata/toLower.bed",
		OutFile:            "testdata/toLower.out.fa",
		IgnoreExtraRegions: true,
		ExpectedFile:       "testdata/expected.toLower.fa"},
}

func TestToLower(t *testing.T) {
	var records []fasta.Fasta
	var regions []Bed
	var file *fileio.EasyWriter
	for _, v := range ToLowerTests {
		records = fasta.Read(v.InFastaFile)
		regions = Read(v.InBedFile)
		ToLower(records, regions, v.IgnoreExtraRegions)
		file = fileio.EasyCreate(v.OutFile)
		fasta.WriteToFileHandle(file, records, 50)
		err := file.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: in ToLower, OutFile did not match ExpectedFile.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var SegregatingSitesTests = []struct {
	InMfaFile       string
	chrom           string
	refStart        int
	OutMfaFile      string
	OutBedFile      string
	ExpectedMfaFile string
	ExpectedBedFile string
}{
	{InMfaFile: "testdata/test.mfa",
		chrom:           "chrTest",
		refStart:        0,
		OutMfaFile:      "testdata/output.fa",
		OutBedFile:      "testdata/output.bed",
		ExpectedMfaFile: "testdata/expected.mfa",
		ExpectedBedFile: "testdata/expected.bed"},
}

func TestSegregatingSites(t *testing.T) {
	var records, answerFa, expectedFa []fasta.Fasta
	var answerBed, expectedBed []Bed
	var mfaFile *fileio.EasyWriter
	var err error
	for _, v := range SegregatingSitesTests {
		records = fasta.Read(v.InMfaFile)
		answerFa, answerBed = SegregatingSites(records, v.chrom, v.refStart)
		mfaFile = fileio.EasyCreate(v.OutMfaFile)
		fasta.WriteToFileHandle(mfaFile, answerFa, 50)
		err = mfaFile.Close()
		exception.PanicOnErr(err)
		Write(v.OutBedFile, answerBed)
		expectedFa = fasta.Read(v.ExpectedMfaFile)
		expectedBed = Read(v.ExpectedBedFile)
		if !fasta.AllAreEqual(answerFa, expectedFa) {
			t.Errorf("Error in mfa part of SegregatingSites.")
		} else {
			err = os.Remove(v.OutMfaFile)
			exception.PanicOnErr(err)
		}
		if !AllAreEqual(answerBed, expectedBed) {
			t.Errorf("Error in bed part of SegregatingSites.")
		} else {
			err = os.Remove(v.OutBedFile)
			exception.PanicOnErr(err)
		}
	}
}
