package main

import (
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

var SimulateSamTests = []struct {
	OutFile                         string
	RefFile                         string
	NumReads                        int
	Coverage                        float64
	ReadLength                      int
	FlatErrorRate                   float64
	InsertLength                    int
	InsertStdDev                    float64
	SetSeed                         int64
	GeometricParam                  float64
	AncientErrorRate                float64
	DeaminationDistribution         string
	ExpectedFile                    string
	ExpectedDeaminationDistribution string
}{
	{OutFile: "testdata/actual.sam",
		RefFile:       "testdata/test.fa",
		NumReads:      100,
		Coverage:      0,
		ReadLength:    150,
		FlatErrorRate: 0,
		InsertLength:  500,
		InsertStdDev:  50,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.sam"},
	{OutFile: "testdata/actual.bam",
		RefFile:       "testdata/test.fa",
		NumReads:      100,
		Coverage:      0,
		ReadLength:    150,
		FlatErrorRate: 0,
		InsertLength:  500,
		InsertStdDev:  50,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.bam"},
	{OutFile: "testdata/10xCoverage.sam",
		RefFile:       "testdata/test.fa",
		NumReads:      100, // this value will be ignored
		Coverage:      10,
		ReadLength:    150,
		FlatErrorRate: 0,
		InsertLength:  50,
		InsertStdDev:  10,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.10xCoverage.sam"},
	{OutFile: "testdata/100xCoverage.sam",
		RefFile:       "testdata/test.fa",
		NumReads:      100, // this value will be ignored
		Coverage:      100,
		ReadLength:    150,
		FlatErrorRate: 0,
		InsertLength:  50,
		InsertStdDev:  10,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.100xCoverage.sam"},
	{OutFile: "testdata/errorTest.LowRate.sam",
		RefFile:       "testdata/errorTest.fa",
		NumReads:      100, // this value will be ignored
		Coverage:      10,
		ReadLength:    50,
		FlatErrorRate: 0.01,
		InsertLength:  100,
		InsertStdDev:  10,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.errorTest.LowRate.sam"},
	{OutFile: "testdata/errorTest.HighRate.sam",
		RefFile:       "testdata/errorTest.fa",
		NumReads:      100, // this value will be ignored
		Coverage:      10,
		ReadLength:    50,
		FlatErrorRate: 0.1,
		InsertLength:  100,
		InsertStdDev:  10,
		SetSeed:       1,
		ExpectedFile:  "testdata/expected.errorTest.HighRate.sam"},
	{OutFile: "testdata/ancientErrorTest.sam",
		RefFile:                         "testdata/test.fa",
		NumReads:                        100, // this value will be ignored
		Coverage:                        10,
		ReadLength:                      50,
		FlatErrorRate:                   0.01,
		InsertLength:                    100,
		InsertStdDev:                    10,
		SetSeed:                         1,
		AncientErrorRate:                0.1,
		GeometricParam:                  0.25,
		DeaminationDistribution:         "testdata/test.deaminationDistribution.txt",
		ExpectedFile:                    "testdata/expected.ancientErrorTest.sam",
		ExpectedDeaminationDistribution: "testdata/expected.deaminationDistribution.txt"},
}

func TestSimulateSam(t *testing.T) {
	var s Settings
	var err error
	var bamA, bamB []sam.Sam
	for _, v := range SimulateSamTests {
		s = Settings{
			OutFile:                 v.OutFile,
			RefFile:                 v.RefFile,
			NumReads:                v.NumReads,
			Coverage:                v.Coverage,
			ReadLength:              v.ReadLength,
			FlatError:               v.FlatErrorRate,
			FragmentLength:          v.InsertLength,
			FragmentStdDev:          v.InsertStdDev,
			SetSeed:                 v.SetSeed,
			AncientErrorRate:        v.AncientErrorRate,
			GeometricParam:          v.GeometricParam,
			DeaminationDistribution: v.DeaminationDistribution,
		}
		simulateSam(s)
		if strings.HasSuffix(s.OutFile, ".bam") {
			bamA, _ = sam.Read(v.OutFile)
			bamB, _ = sam.Read(v.ExpectedFile)
			if len(bamA) != len(bamB) {
				t.Error("Error: problem simulating bam")
			} else {
				for i := range bamA {
					if !sam.Equal(bamA[i], bamB[i]) {
						t.Error("Error: problem simulating bam")
						fmt.Printf("Observed:\n%s\n", sam.ToString(bamA[i]))
						fmt.Printf("Expected:\n%s\n", sam.ToString(bamB[i]))
					}
				}
			}
			if !t.Failed() {
				err = os.Remove(v.OutFile)
				exception.PanicOnErr(err)
			}
		} else {
			if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
				t.Errorf("Error: Problem simulating sam.")
			} else {
				err = os.Remove(v.OutFile)
				exception.PanicOnErr(err)
			}
		}

		if v.ExpectedDeaminationDistribution != "" {
			if !fileio.AreEqual(v.DeaminationDistribution, v.ExpectedDeaminationDistribution) {
				t.Errorf("Error: deamination distribution file did not match expected.")
			}
		}

	}
}
