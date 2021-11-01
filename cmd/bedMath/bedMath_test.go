package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var BedMathTests = []struct {
	Afile        string
	Bfile        string
	Outfile      string
	ExpectedFile string
	Op           string
}{
	{"testdata/testA.bed", "testdata/testB.bed", "testdata/test.Add.bed", "testdata/expected.Add.bed", "Add"},
	{"testdata/testA.bed", "testdata/testB.bed", "testdata/test.Sub.bed", "testdata/expected.Sub.bed", "Subtract"},
	{"testdata/testA.bed", "testdata/testB.bed", "testdata/test.Mult.bed", "testdata/expected.Mult.bed", "Multiply"},
	{"testdata/testA.bed", "testdata/testB.bed", "testdata/test.Divide.bed", "testdata/expected.Divide.bed", "Divide"},
}

func TestBedMath(t *testing.T) {
	var err error
	var s Settings
	for _, v := range BedMathTests {
		s = Settings{
			Afile:   v.Afile,
			Bfile:   v.Bfile,
			OutFile: v.Outfile,
			Op:      v.Op,
		}
		bedMath(s)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in bedMath.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}
