package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var DunnIndexTests = []struct {
	InputBed     string
	InputFa      string
	InputGroup   string
	OutFile      string
	ExpectedFile string
	Realign      bool
}{
	{"testdata/test.bed", "testdata/test.fa", "testdata/groups.list", "testdata/tmp.bed", "testdata/expected.bed", false},
	{"testdata/test.realign.bed", "testdata/test.realign.fa", "testdata/groups.list", "testdata/tmp.realign.bed", "testdata/expected.realign.bed", true},
}

func TestDunnIndex(t *testing.T) {
	var err error
	for _, v := range DunnIndexTests {
		dunnIndex(v.InputBed, v.InputFa, v.InputGroup, v.Realign, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in dunnIndex. Output did not match expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
