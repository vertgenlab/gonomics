package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var BranchLengthsMultiFaBedTests = []struct {
	Chrom                      string
	InFaFile                   string
	InBedFile                  string
	VelLengthBedFile           string
	InitialLengthBedFile       string
	NumUngappedSitesBedFile    string
	QOutFile                   string
	SearchSpaceBed             string
	SearchSpaceProportion      float64
	UseSnpDistance             bool
	Verbose                    bool
	Epsilon                    float64
	AllowNegative              bool
	ZeroDistanceWeightConstant float64
	VelLengthBedExpected       string
	InitialLengthBedExpected   string
	NumUngappedSitesExpected   string
	QoutExpected               string
	CavalliSforzaQ             bool
}{
	{"chr1",
		"testdata/test.fa",
		"testdata/test.in.bed",
		"testdata/test.Vel.bed",
		"testdata/test.Initial.bed",
		"testdata/test.UngappedSites.bed",
		"testdata/qOut.bed",
		"",
		0.5,
		false,
		false,
		1e-8,
		false,
		1000,
		"testdata/expected.Vel.bed",
		"testdata/expected.Initial.bed",
		"testdata/expected.NumUngapped.bed",
		"testdata/QoutExpected.bed",
		false,
	},
}

func TestBranchLengthsMultiFaBed(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BranchLengthsMultiFaBedTests {
		s = Settings{
			Chrom:                      v.Chrom,
			InFaFile:                   v.InFaFile,
			InBedFile:                  v.InBedFile,
			VelLengthBedFile:           v.VelLengthBedFile,
			InitialLengthBedFile:       v.InitialLengthBedFile,
			NumUngappedSitesBedFile:    v.NumUngappedSitesBedFile,
			QOutFile:                   v.QOutFile,
			SearchSpaceBed:             v.SearchSpaceBed,
			SearchSpaceProportion:      v.SearchSpaceProportion,
			UseSnpDistance:             v.UseSnpDistance,
			Verbose:                    v.Verbose,
			Epsilon:                    v.Epsilon,
			AllowNegative:              v.AllowNegative,
			ZeroDistanceWeightConstant: v.ZeroDistanceWeightConstant,
			CavalliSforzaEdwardsQ:      v.CavalliSforzaQ,
		}
		branchLengthsMultiFaBed(s)
		if !fileio.AreEqual(v.VelLengthBedFile, v.VelLengthBedExpected) {
			t.Errorf("Error in branchLengthsMultiFaBed. VelLengthBedFile did not match expected.")
		} else {
			err = os.Remove(v.VelLengthBedFile)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.InitialLengthBedFile, v.InitialLengthBedExpected) {
			t.Errorf("Error in branchLengthsMultiFaBed. InitialLengthBedFile did not match expected.")
		} else {
			err = os.Remove(v.InitialLengthBedFile)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.NumUngappedSitesBedFile, v.NumUngappedSitesExpected) {
			t.Errorf("Error in branchLengthsMultiFaBed. NumUngappedSitesBedFile did not match expected.")
		} else {
			err = os.Remove(v.NumUngappedSitesBedFile)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.QOutFile, v.QoutExpected) {
			t.Errorf("Error in branchLengthsMultiFaBed. QoutFile did not match expected.")
		} else {
			err = os.Remove(v.QOutFile)
			exception.PanicOnErr(err)
		}
	}
}
