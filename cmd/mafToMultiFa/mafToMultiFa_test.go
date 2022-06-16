package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var mafToMultiFaTests = []struct {
	mafFile                string
	refFile                string
	speciesFile            string
	outMultiFaFileExpected string
	noMask                 bool
}{
	{"testdata/test1.maf", "testdata/test.ref.fa", "testdata/test.species.list", "testdata/test.out.fa", false},
	{"testdata/test2.maf", "testdata/test.ref.fa", "testdata/test.species.list", "testdata/test.out.fa", true},
}

func TestMafToMultiFa(t *testing.T) {
	var err error
	for _, v := range mafToMultiFaTests {
		mafToMultiFa(v.mafFile, v.refFile, v.speciesFile, "outFile_tmp.fa", v.noMask)

		if !fileio.AreEqual("outFile_tmp.fa", v.outMultiFaFileExpected) {
			t.Errorf("Error in mafToFa")
		} else {
			err = os.Remove("outFile_tmp.fa")
		}
		exception.PanicOnErr(err)
	}
}
