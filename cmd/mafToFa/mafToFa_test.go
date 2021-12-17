package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var mafToFaTests = []struct {
	mafFile	string
	refFile	string
	speciesFile string
	outFile_expected string
	noMask bool
}{
	{"testdata/test1.maf", "testdata/test.ref.fa", "testdata/test.species.list", "testdata/test.out.fa", false},
	{"testdata/test2.maf", "testdata/test.ref.fa", "testdata/test.species.list", "testdata/test.out.fa", true},
}

func TestMafToFa(t *testing.T) {
	var err error
	for _, v := range mafToFaTests {
		mafToFa(v.mafFile, v.refFile, v.speciesFile, "outFile_tmp.mfa", v.noMask)

		if !fileio.AreEqual("outFile_tmp.mfa", v.outFile_expected) {
			t.Errorf("Error in mafToFa")
		} else {
			err = os.Remove("outFile_tmp.mfa")
		}
		exception.PanicOnErr(err)
	}
}
