package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var SimulateEvolTests = []struct {
	InputFastaFile string
	OutFile string
	TreeFile string
	SimOutFile string
	GenePredFile string
	BranchLength float64
	Lambda float64
	PropIndel float64
	GcContent float64
	SetSeed int64
	VcfOutFile string
	ExpectedOutFile string
	ExpectedVcfFile string
}{
	{InputFastaFile: "testdata/rand.fa",
		OutFile: "testdata/tmp.out.fa",
	TreeFile: "",
	SimOutFile: "",
	GenePredFile: "",
	BranchLength: 0.1,
	Lambda: 1,
	PropIndel: 0.2,
	GcContent: 0.42,
	SetSeed: -1,
	VcfOutFile: "testdata/tmp.vcf",
	ExpectedVcfFile: "testdata/expected.branchLength.vcf",
	ExpectedOutFile: "testdata/expected.branchLength.fa"},
}

func TestSimulateEvol(t *testing.T) {
	var err error
	var s Settings
	for _, v := range SimulateEvolTests {
		s = Settings {
			FastaFile: v.InputFastaFile,
			TreeFile: v.TreeFile,
			LeafOutFile: v.OutFile,
			SimOutFile: v.SimOutFile,
			GenePredFile: v.GenePredFile,
			BranchLength: v.BranchLength,
			Lambda: v.Lambda,
			PropIndel: v.PropIndel,
			GcContent: v.GcContent,
			SetSeed: v.SetSeed,
			VcfOutFile: v.VcfOutFile,
		}
		SimulateEvol(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedOutFile) {
			t.Errorf("Error in SimulateEvol. OutFile was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if v.VcfOutFile != "" {
			if !fileio.AreEqual(v.VcfOutFile, v.ExpectedVcfFile) {
				t.Errorf("Error in SimulateEvol. Output Vcf file was not as expected.")
			} else {
				err = os.Remove(v.VcfOutFile)
				exception.PanicOnErr(err)
			}
		}
	}
}