package net

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestReadAndWrite(t *testing.T) {
	var outfile, infile string = "testdata/out.net", "testdata/test.in.net"
	in, mp := Read(infile)
	Write(outfile, in, mp)
	if !fileio.AreEqual(outfile, infile) {
		t.Errorf("Error in Net reading and writing\n")
	}
	err := os.Remove(outfile)
	exception.PanicOnErr(err)
}

func TestGoReadToChan(t *testing.T) {
	var n []Net
	var outfile, infile string = "testdata/out.net", "testdata/test.in.net"
	inChan, mp := GoReadToChan(infile)
	for i := range inChan {
		n = append(n, i)
	}
	Write(outfile, n, mp)
	if !fileio.AreEqual(outfile, infile) {
		t.Errorf("Error in GoReadToChan\n")
	} else {
		err := os.Remove(outfile)
		exception.PanicOnErr(err)
	}
}
