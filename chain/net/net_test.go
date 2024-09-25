package net

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestReadAndWrite(t *testing.T) {
	var out, exp string = "testdata/out.net", "testdata/exp.net"
	in, mp := Read("testdata/test.in.net")
	Write(out, in, mp)
	if !fileio.AreEqual(out, exp) {
		t.Errorf("Error in Net reading and writing\n")
	}
	err := os.Remove(out)
	exception.PanicOnErr(err)
}
