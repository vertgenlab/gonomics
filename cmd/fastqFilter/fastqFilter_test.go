package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
)

var FastqFilterTests = []struct {
	inputFile        string
	outputFile       string
	expectedFile     string
	R1InFile         string
	R2InFile         string
	R1OutFile        string
	R2OutFile        string
	R1ExpectedFile   string
	R2ExpectedFile   string
	PairedEnd        bool
	SubSet           float64
	SetSeed          int64
	MinSize          int
	MaxSize          int
	RetainNamesList  string
	DiscardNamesList string
	CollapseUmi      bool
	BarcodeLength    int
	UmiLength        int
}{
	{"../../fastq/testdata/test.fastq", "tmpOut.fastq", "testdata/expectedReadWrite.fastq", "", "", "", "", "", "", false, 1, 10, 0, numbers.MaxInt, "", "", false, 16, 12},
	{"../../fastq/testdata/test.fastq", "tmpOut.fastq", "testdata/expectedHalf.fastq", "", "", "", "", "", "", false, 0.5, 10, 0, numbers.MaxInt, "", "", false, 16, 12},
	{"", "", "", "../../fastq/testdata/simReads_R1.fq", "../../fastq/testdata/simReads_R2.fq", "tmpR1.fastq", "tmpR2.fastq", "testdata/expectedR1ReadWrite.fastq", "testdata/expectedR2ReadWrite.fastq", true, 1, 10, 0, numbers.MaxInt, "", "", false, 16, 12}, //~/go/bin/fastqFilter -pairedEnd -setSeed 10 ../../fastq/testdata/simReads_R1.fq ../../fastq/testdata/simReads_R2.fq testdata/expectedR1ReadWrite.fastq testdata/expectedR2ReadWrite.fastq
	{"", "", "", "../../fastq/testdata/simReads_R1.fq", "../../fastq/testdata/simReads_R2.fq", "tmpR1.fastq", "tmpR2.fastq", "testdata/expectedR1Half.fastq", "testdata/expectedR2Half.fastq", true, 0.5, 10, 0, numbers.MaxInt, "", "", false, 16, 12},         //~/go/bin/fastqFilter -pairedEnd -setSeed 10 -subSet 0.5 ../../fastq/testdata/simReads_R1.fq ../../fastq/testdata/simReads_R2.fq testdata/expectedR1Half.fastq testdata/expectedR2Half.fastq
	{"", "", "", "testdata/UmiTest_R1.fastq", "testdata/UmiTest_R2.fastq", "tmpR1.fastq", "tmpR2.fastq", "testdata/expectedUmi_R1.fastq", "testdata/expectedUmi_R2.fastq", true, 1, 10, 0, numbers.MaxInt, "", "", true, 16, 12},
	{"../../fastq/testdata/test.fastq", "tmpOut.fastq", "testdata/expectedNamesFilter.fastq", "", "", "", "", "", "", false, 1, 10, 0, numbers.MaxInt, "testdata/namesList.txt", "", false, 16, 12},        //~/go/bin/fastqFilter -setSeed 10 -retainNamesList testdata/namesList.txt ../../fastq/testdata/test.fastq testdata/expectedNamesFilter.fastq
	{"../../fastq/testdata/test.fastq", "tmpOut.fastq", "testdata/expectedDiscardNamesFilter.fastq", "", "", "", "", "", "", false, 1, 10, 0, numbers.MaxInt, "", "testdata/namesList.txt", false, 16, 12}, //~/go/bin/fastqFilter -setSeed 10 -discardNamesList testdata/namesList.txt ../../fastq/testdata/test.fastq testdata/expectedDiscardNamesFilter.fastq
}

func TestFastqFilter(t *testing.T) {
	var err error
	var s Settings
	for _, v := range FastqFilterTests {
		s = Settings{
			InFile:           v.inputFile,
			OutFile:          v.outputFile,
			R1InFile:         v.R1InFile,
			R2InFile:         v.R2InFile,
			R1OutFile:        v.R1OutFile,
			R2OutFile:        v.R2OutFile,
			PairedEnd:        v.PairedEnd,
			SubSet:           v.SubSet,
			SetSeed:          v.SetSeed,
			MinSize:          v.MinSize,
			MaxSize:          v.MaxSize,
			RetainNamesList:  v.RetainNamesList,
			DiscardNamesList: v.DiscardNamesList,
			CollapseUmi:      v.CollapseUmi,
			BarcodeLength:    v.BarcodeLength,
			UmiLength:        v.UmiLength,
		}
		fastqFilter(s)
		if v.PairedEnd {
			if !fileio.AreEqual(v.R1OutFile, v.R1ExpectedFile) {
				t.Errorf("Error in fastqFilter, paired reads, read 1.")
			}
			if !fileio.AreEqual(v.R2OutFile, v.R2ExpectedFile) {
				t.Errorf("Error in fastqFilter, paired reads, read 2.")
			}
			err = os.Remove(v.R1OutFile)
			if err != nil {
				common.ExitIfError(err)
			}
			err = os.Remove(v.R2OutFile)
			if err != nil {
				common.ExitIfError(err)
			}
		} else {
			if !fileio.AreEqual(v.outputFile, v.expectedFile) {
				t.Errorf("Error in fastqFilter, unpaired reads.")
			}
			err = os.Remove(v.outputFile)
			if err != nil {
				common.ExitIfError(err)
			}
		}
	}
}
