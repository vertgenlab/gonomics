package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"os"
	"strings"
	"testing"
)

var SimulateSamTests = []struct {
	OutFile      string
	RefFile      string
	NumReads     int
	ReadLength   int
	InsertLength int
	InsertStdDev float64
	SetSeed      int64
	ExpectedFile string
}{
	{OutFile: "testdata/actual.sam",
		RefFile:      "testdata/test.fa",
		NumReads:     100,
		ReadLength:   150,
		InsertLength: 500,
		InsertStdDev: 50,
		SetSeed:      1,
		ExpectedFile: "testdata/expected.sam"},
	{OutFile: "testdata/actual.bam",
		RefFile:      "testdata/test.fa",
		NumReads:     100,
		ReadLength:   150,
		InsertLength: 500,
		InsertStdDev: 50,
		SetSeed:      1,
		ExpectedFile: "testdata/expected.bam"},
}

func TestSimulateSam(t *testing.T) {
	var s Settings
	var err error
	var bamA, bamB []sam.Sam
	for _, v := range SimulateSamTests {
		s = Settings{
			OutFile:        v.OutFile,
			RefFile:        v.RefFile,
			NumReads:       v.NumReads,
			ReadLength:     v.ReadLength,
			FragmentLength: v.InsertLength,
			FragmentStdDev: v.InsertStdDev,
			SetSeed:        v.SetSeed,
		}
		simulateSam(s)
		if strings.HasSuffix(s.OutFile, ".bam") {
			bamA, _ = sam.Read(v.OutFile)
			bamB, _ = sam.Read(v.ExpectedFile)
			if len(bamA) != len(bamB) {
				t.Error("problem simulating bam")
			} else {
				for i := range bamA {
					if !sam.Equal(bamA[i], bamB[i]) {
						t.Error("problem simulating bam")
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
				t.Errorf("Problem simulating sam.")
			} else {
				err = os.Remove(v.OutFile)
				exception.PanicOnErr(err)
			}
		}
	}
}
