package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var DunnIndexTests = []struct {
	InputBed     string
	InputFa      string
	InputGroup   string
	OutFile      string
	ExpectedFile string
}{
	{"testdata/test.bed", "testdata/test.fa", "testdata/groups.list", "testdata/tmp.bed", "testdata/expected.bed"},
}

func TestDunnIndex(t *testing.T) {
	var err error
	for _, v := range DunnIndexTests {
		dunnIndex(v.InputBed, v.InputFa, v.InputGroup, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in dunnIndex. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
