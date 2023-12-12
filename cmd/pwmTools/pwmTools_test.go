package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var filterTests = []struct {
	InFile       string
	OutFile      string
	MatrixType   string
	MaxLength    int
	MinLength    int
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/jaspar.vertebrate.filterMin10Max15.txt",
		MatrixType:   "Frequency",
		MinLength:    10,
		MaxLength:    15,
		ExpectedFile: "testdata/expected.filterMin10Max15.txt",
	},
}

func TestPwmFilter(t *testing.T) {
	var err error
	var s FilterSettings
	for _, v := range filterTests {
		s = FilterSettings{
			InFile:     v.InFile,
			OutFile:    v.OutFile,
			MatrixType: v.MatrixType,
			MaxLength:  v.MaxLength,
			MinLength:  v.MinLength,
		}
		pwmFilter(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmFilter outFile is not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var InfoTests = []struct {
	InFile       string
	OutFile      string
	MatrixType   string
	PseudoCounts float64
	GcContent    float64
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/jaspar.vertebrate.info.txt",
		MatrixType:   "Frequency",
		PseudoCounts: 0.1,
		GcContent:    0.5,
		ExpectedFile: "testdata/expected.info.txt",
	},
}

func TestPwmInfo(t *testing.T) {
	var err error
	var s InfoSettings
	for _, v := range InfoTests {
		s = InfoSettings{
			InFile:       v.InFile,
			OutFile:      v.OutFile,
			MatrixType:   v.MatrixType,
			PseudoCounts: v.PseudoCounts,
			GcContent:    v.GcContent,
		}
		pwmInfo(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmInfo outFile is not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
