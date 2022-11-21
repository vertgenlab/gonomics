package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultiFaToAxtTests = []struct {
	InFile        string
	RName         string
	QName         string
	OutFile       string
	ExpecctedFile string
}{
	{"testdata/test.fa", "chr1", "chr1", "testdata/tmp.axt", "testdata/expected.axt"},
}

func TestMultiFaToAxt(t *testing.T) {
	var err error
	for _, v := range MultiFaToAxtTests {
		multiFaToAxt(v.InFile, v.RName, v.QName, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpecctedFile) {
			t.Errorf("Error in multiFatoAxt. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
