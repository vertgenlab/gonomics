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

var FormatTests = []struct {
	InFile       string
	OutFile      string
	InType       string
	OutType      string
	PseudoCount  float64
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.jaspar.ppm.txt",
		InType:       "Frequency",
		OutType:      "Probability",
		PseudoCount:  0,
		ExpectedFile: "testdata/expected.jaspar.ppm.txt",
	},
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.highPseudo.jaspar.ppm.txt",
		InType:       "Frequency",
		OutType:      "Probability",
		PseudoCount:  40,
		ExpectedFile: "testdata/expected.highPseudo.jaspar.ppm.txt",
	},
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.jaspar.pwm.txt",
		InType:       "Frequency",
		OutType:      "Weight",
		PseudoCount:  0.2,
		ExpectedFile: "testdata/expected.jaspar.pwm.txt",
	},
}

func TestPwmFormat(t *testing.T) {
	var err error
	var s FormatSettings
	for _, v := range FormatTests {
		s = FormatSettings{
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			InType:      v.InType,
			OutType:     v.OutType,
			PseudoCount: v.PseudoCount,
		}
		pwmFormat(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmFormat outFile is not as expected.")
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
	Pseudocounts float64
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/jaspar.vertebrate.info.txt",
		MatrixType:   "Frequency",
		Pseudocounts: 0.1,
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
			Pseudocounts: v.Pseudocounts,
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
