package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultiFaToChainTests = []struct {
	InFile string
	tName string
	qName string
	OutFile string
	ExpectedFile string
}{
	{InFile: "testdata/test.fa",
		tName: "chr22",
	qName: "chr22",
	OutFile: "testdata/tmp.chain",
	ExpectedFile: "testdata/expected.chain",},
}

func TestMultiFaToChain(t *testing.T) {
	var err error
	for _, v := range MultiFaToChainTests {
		multiFaToChain(v.InFile, v.tName, v.qName, v.OutFile)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in MultiFaToChain. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

