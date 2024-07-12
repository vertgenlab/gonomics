package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
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
			fileio.EasyRemove(v.OutFile)
		}
	}
}

var FormatTests = []struct {
	InFile       string
	OutFile      string
	InType       string
	OutType      string
	PseudoCount  float64
	GcContent    float64
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.jaspar.ppm.txt",
		InType:       "Frequency",
		OutType:      "Probability",
		PseudoCount:  0,
		GcContent:    0.5,
		ExpectedFile: "testdata/expected.jaspar.ppm.txt",
	},
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.highPseudo.jaspar.ppm.txt",
		InType:       "Frequency",
		OutType:      "Probability",
		PseudoCount:  40,
		GcContent:    0.5,
		ExpectedFile: "testdata/expected.highPseudo.jaspar.ppm.txt",
	},
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/tmp.jaspar.pwm.txt",
		InType:       "Frequency",
		OutType:      "Weight",
		PseudoCount:  0.2,
		GcContent:    0.5,
		ExpectedFile: "testdata/expected.jaspar.pwm.txt",
	},
}

func TestPwmFormat(t *testing.T) {
	var s FormatSettings
	for _, v := range FormatTests {
		s = FormatSettings{
			InFile:      v.InFile,
			OutFile:     v.OutFile,
			InType:      v.InType,
			OutType:     v.OutType,
			PseudoCount: v.PseudoCount,
			GcContent:   v.GcContent,
		}
		pwmFormat(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmFormat outFile is not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}

var InfoTests = []struct {
	InFile       string
	OutFile      string
	MatrixType   string
	PseudoCounts float64
	GcContent    float64
	Threshold    float64
	ExpectedFile string
}{
	{InFile: "testdata/jaspar.vertebrate.txt.gz",
		OutFile:      "testdata/jaspar.vertebrate.info.txt",
		MatrixType:   "Frequency",
		PseudoCounts: 0.1,
		GcContent:    0.5,
		Threshold:    0.8,
		ExpectedFile: "testdata/expected.info.txt",
	},
}

func TestPwmInfo(t *testing.T) {
	var s InfoSettings
	for _, v := range InfoTests {
		s = InfoSettings{
			InFile:       v.InFile,
			OutFile:      v.OutFile,
			MatrixType:   v.MatrixType,
			PseudoCounts: v.PseudoCounts,
			GcContent:    v.GcContent,
			Threshold:    v.Threshold,
		}
		pwmInfo(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmInfo outFile is not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}

var PwmShuffleTests = []struct {
	InFile       string
	OutFile      string
	NumShuffle   int
	SetSeed      int64
	ExpectedFile string
}{
	{InFile: "testdata/firstSix.jaspar.pwm.txt",
		OutFile:      "testdata/test.firstSix.shuffle.pwm.txt",
		NumShuffle:   10,
		SetSeed:      13,
		ExpectedFile: "testdata/expected.firstSix.shuffle.pwm.txt",
	},
}

func TestPwmShuffle(t *testing.T) {
	var s ShuffleSettings
	for _, v := range PwmShuffleTests {
		s = ShuffleSettings{
			InFile:     v.InFile,
			OutFile:    v.OutFile,
			NumShuffle: v.NumShuffle,
			SetSeed:    v.SetSeed,
		}
		pwmShuffle(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: pwmShuffle outFile is not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
	}
}
