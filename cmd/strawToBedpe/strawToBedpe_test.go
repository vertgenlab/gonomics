package main

import (
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var StrawToBedPeTests = []struct {
	FileList                          string
	OutFile                           string
	BinSize                           int
	RStart                            float64
	PStart                            float64
	RStep                             float64
	PStep                             float64
	MinCutoff                         int
	MinBinDistance                    int
	Fdr                               float64
	FitStatsFile                      string
	ContactScoreFile                  string
	Expected                          string
	ExpectedFitStats                  string
	ExpectedContactScoreFile          string
	MaxContactScoreInDistributionFile int
	MaxBinDistanceInDistributionFile  int
}{
	{FileList: "testdata/fileList.txt",
		OutFile:                           "testdata/out.bedpe",
		BinSize:                           5000,
		RStart:                            1.0,
		PStart:                            0.5,
		RStep:                             0.001,
		PStep:                             0.001,
		MinBinDistance:                    2,
		MinCutoff:                         10,
		ContactScoreFile:                  "testdata/out.contactScoreFile.txt",
		FitStatsFile:                      "testdata/out.FitStats.txt",
		Fdr:                               0.05,
		Expected:                          "testdata/expected.out.bedpe",
		ExpectedFitStats:                  "testdata/expected.FitStats.txt",
		ExpectedContactScoreFile:          "testdata/expected.contactScoreFile.txt",
		MaxContactScoreInDistributionFile: 100,
		MaxBinDistanceInDistributionFile:  -1,
	},
	{FileList: "testdata/fileList.txt",
		OutFile:          "testdata/out.lowCutoff.bedpe",
		BinSize:          5000,
		RStart:           1.0,
		PStart:           0.5,
		RStep:            0.001,
		PStep:            0.001,
		MinCutoff:        2,
		FitStatsFile:     "testdata/out.FitStats.lowCoverage.txt",
		Fdr:              0.05,
		Expected:         "testdata/expected.out.lowCutoff.bedpe",
		ExpectedFitStats: "testdata/expected.FitStats.lowCoverage.txt",
	},
}

func TestStrawToBedpe(t *testing.T) {
	var err error
	var s Settings
	for _, v := range StrawToBedPeTests {
		s = Settings{
			FileList:                          v.FileList,
			OutFile:                           v.OutFile,
			BinSize:                           v.BinSize,
			RStart:                            v.RStart,
			PStart:                            v.PStart,
			RStep:                             v.RStep,
			PStep:                             v.PStep,
			MinCutoff:                         v.MinCutoff,
			MinBinDistance:                    v.MinBinDistance,
			Fdr:                               v.Fdr,
			FitStatsFile:                      v.FitStatsFile,
			ContactScoreFile:                  v.ContactScoreFile,
			MaxBinDistanceInDistributionFile:  v.MaxBinDistanceInDistributionFile,
			MaxContactScoreInDistributionFile: v.MaxContactScoreInDistributionFile,
		}
		strawToBedpe(s)
		if !bedpe.AllAreEqual(bedpe.Read(v.Expected), bedpe.Read(v.OutFile)) {
			t.Errorf("outFile: %s did not match expected file: %s.", v.OutFile, v.Expected)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
		if v.ExpectedFitStats != "" && !fileio.AreEqual(v.ExpectedFitStats, v.FitStatsFile) {
			t.Errorf("Error: Fit stats file did not match expected.\n")
		} else if v.ExpectedFitStats != "" {
			err = os.Remove(v.FitStatsFile)
			exception.PanicOnErr(err)
		}
		if v.ExpectedContactScoreFile != "" && !fileio.AreEqual(v.ExpectedContactScoreFile, v.ContactScoreFile) {
			t.Errorf("Error: ContactScore file did not match expected.\n")
		} else if v.ExpectedContactScoreFile != "" {
			err = os.Remove(v.ContactScoreFile)
			exception.PanicOnErr(err)
		}
	}
}
