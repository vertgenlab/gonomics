package net

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestToBed(t *testing.T) {
	exp := "testdata/exp.NTB.bed"
	out := "testdata/out.NTB.bed"
	nets, _ := Read("testdata/test.in.net")
	beds := ToBed(nets)
	bed.Write(out, beds)
	if !fileio.AreEqual(exp, out) {
		t.Errorf("Error in ToBed")
	} else {
		err := os.Remove(out)
		exception.PanicOnErr(err)
	}
}
