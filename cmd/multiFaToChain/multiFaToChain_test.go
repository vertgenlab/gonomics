package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var MultiFaToChainTests = []struct {
	InFile       string
	tName        string
	qName        string
	OutFile      string
	SwapTandQ    bool
	ExpectedFile string
}{
	{InFile: "testdata/test.fa",
		tName:        "chr22",
		qName:        "chr22",
		OutFile:      "testdata/tmp.chain",
		SwapTandQ:    false,
		ExpectedFile: "testdata/expected.chain"},
	{InFile: "testdata/test.fa",
		tName:        "chr22",
		qName:        "chr22",
		OutFile:      "testdata/tmp.swap.chain",
		SwapTandQ:    true,
		ExpectedFile: "testdata/expected.swap.chain"},
}

func TestMultiFaToChain(t *testing.T) {
	var err error
	for _, v := range MultiFaToChainTests {
		multiFaToChain(v.InFile, v.tName, v.qName, v.OutFile, v.SwapTandQ)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in MultiFaToChain. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
