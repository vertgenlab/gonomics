package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ProximityBlockVcfTests = []struct {
	InFile string
	OutFile string
	ExpectedFile string
	Distance int
	SetSeed int64
}{
	{"testdata/test.vcf", "testdata/tmp.vcf", "testdata/expectedSeedMinus1.vcf", 10, -1},
	{"testdata/test.vcf", "testdata/tmp.vcf", "testdata/expectedSeed10.vcf", 10, 10},
}

func TestProximityBlockVcf(t *testing.T) {
	var err error
	for _, v := range ProximityBlockVcfTests {
		proximityBlockVcf(v.InFile, v.OutFile, v.Distance, false, v.SetSeed)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in proximityBlockVcf. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}