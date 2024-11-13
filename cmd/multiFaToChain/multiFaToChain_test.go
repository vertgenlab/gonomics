package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var MultiFaToChainTests = []struct {
	InFile       string
	tName        string
	qName        string
	OutFile      string
	SwapTandQ    bool
	ExpectedFile string
	QuerySeqName string
}{
	{InFile: "testdata/test.fa",
		tName:        "chr22",
		qName:        "chr22",
		OutFile:      "testdata/tmp.chain",
		SwapTandQ:    false,
		ExpectedFile: "testdata/expected.chain",
		QuerySeqName: "",
	},
	{InFile: "testdata/test.fa",
		tName:        "chr22",
		qName:        "chr22",
		OutFile:      "testdata/tmp.swap.chain",
		SwapTandQ:    true,
		ExpectedFile: "testdata/expected.swap.chain",
		QuerySeqName: "",
	},
	{InFile: "testdata/test.ThreeWay.fa",
		tName:        "chr22",
		qName:        "chr22",
		OutFile:      "testdata/tmp.threeWay.chain",
		SwapTandQ:    false,
		ExpectedFile: "testdata/expected.chain",
		QuerySeqName: "hca",
	},
}

func TestMultiFaToChain(t *testing.T) {
	var err error
	var s Settings
	for _, v := range MultiFaToChainTests {
		s = Settings{
			InFile:       v.InFile,
			TName:        v.tName,
			QName:        v.qName,
			OutFile:      v.OutFile,
			SwapTandQ:    v.SwapTandQ,
			QuerySeqName: v.QuerySeqName,
		}
		multiFaToChain(s)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in MultiFaToChain. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
