package main

import (
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var WithIndelsTests = []struct {
	InputFastaFile  string
	OutFile         string
	BranchLength    float64
	Lambda          float64
	PropIndels      float64
	GcContent       float64
	SetSeed         int64
	TransitionBias  float64
	VcfOutFile      string
	ExpectedOutFile string
	ExpectedVcfFile string
	QName           string
}{
	{InputFastaFile: "testdata/rand.fa.gz",
		OutFile:         "testdata/tmp.out.fa",
		BranchLength:    0.1,
		Lambda:          1,
		PropIndels:      0.2,
		GcContent:       0.42,
		SetSeed:         -1,
		TransitionBias:  1,
		VcfOutFile:      "testdata/tmp.vcf",
		ExpectedVcfFile: "testdata/expected.branchLength.vcf.gz",
		ExpectedOutFile: "testdata/expected.branchLength.fa.gz",
		QName:           "sim"},
}

func TestSimulateEvol(t *testing.T) {

	var s WithIndelsSettings
	for _, v := range WithIndelsTests {
		s = WithIndelsSettings{
			FastaFile:      v.InputFastaFile,
			OutFile:        v.OutFile,
			BranchLength:   v.BranchLength,
			Lambda:         v.Lambda,
			PropIndels:     v.PropIndels,
			GcContent:      v.GcContent,
			SetSeed:        v.SetSeed,
			VcfOutFile:     v.VcfOutFile,
			TransitionBias: v.TransitionBias,
			QName:          v.QName,
		}
		WithIndels(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedOutFile) {
			t.Errorf("Error in SimulateEvol withIndels. OutFile was not as expected.")
		} else {
			fileio.EasyRemove(v.OutFile)
		}
		if v.VcfOutFile != "" {

			if !fileio.AreEqual(v.VcfOutFile, v.ExpectedVcfFile) {
				t.Errorf("Error in SimulateEvol withIndels. Output Vcf file was not as expected.")
			} else {
				fileio.EasyRemove(v.VcfOutFile)
			}
		}
	}
}
