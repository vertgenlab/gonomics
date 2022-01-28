package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var FaFindFastTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	WindowSize   int
	ChromName    string
	RemoveN      bool
}{
	{"testdata/test_indel.fa", "testdata/tmp.out.bed", "testdata/expected.bed", 10, "chr1", false},
	{"testdata/test_indel.fa", "testdata/tmp.noN.bed", "testdata/expected.noN.bed", 10, "chr1", true},
}

func TestFaFindFast(t *testing.T) {
	var err error
	for _, v := range FaFindFastTests {
		faFindFast(v.InFile, v.OutFile, v.WindowSize, v.ChromName, v.RemoveN, false)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in faFindFast. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
