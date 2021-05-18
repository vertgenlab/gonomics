package main

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FastqFormatTests = []struct {
	InFile            string
	OutFile           string
	ExpectedOutFile   string
	R1InFile          string
	R2InFile          string
	R1OutFile         string
	R2OutFile         string
	ExpectedR1OutFile string
	ExpectedR2OutFile string
	PairedEnd         bool
	SingleCell        bool
	BarcodeLength     int
	UmiLength         int
}{
	{"", "", "", "testdata/TestR1.fastq", "testdata/TestR2.fastq", "tmpR1.fastq", "tmpR2.fastq", "testdata/ExpectedR1.fastq", "testdata/ExpectedR2.fastq", true, true, 16, 12},
}

func TestFastqFormat(t *testing.T) {
	var err error
	var s Settings
	for _, v := range FastqFormatTests {
		s = Settings{
			InFile:        s.InFile,
			OutFile:       s.OutFile,
			R1InFile:      v.R1InFile,
			R2InFile:      v.R2InFile,
			R1OutFile:     v.R1OutFile,
			R2OutFile:     v.R2OutFile,
			PairedEnd:     v.PairedEnd,
			SingleCell:    v.SingleCell,
			BarcodeLength: v.BarcodeLength,
			UmiLength:     v.UmiLength,
		}
		fastqFormat(s)
		if v.PairedEnd {
			if !fileio.AreEqual(v.R1OutFile, v.ExpectedR1OutFile) {
				t.Errorf("Error in fastqFormat, paired reads, read 1.")
			}
			if !fileio.AreEqual(v.R2OutFile, v.ExpectedR2OutFile) {
				t.Errorf("Error in fastqFormat, paired reads, read 2.")
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
			if !fileio.AreEqual(v.OutFile, v.ExpectedOutFile) {
				t.Errorf("Error in fastqFormat, unpaired reads.")
			}
			err = os.Remove(v.OutFile)
			if err != nil {
				common.ExitIfError(err)
			}
		}
	}
}
