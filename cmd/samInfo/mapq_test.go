package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestMapQ(t *testing.T) {
	exp := "testdata/mapq/exp.hist.txt"
	obs := "testdata/mapq/out.hist.txt"
	s := mapqSettings{
		InFile:  "testdata/mapq/test1.sam",
		OutFile: obs,
	}

	mapq(s)
	if !fileio.AreEqual(exp, obs) {
		t.Errorf("Error mapq: The exptected and observed files are not equal\n")
	} else {
		exception.PanicOnErr(os.Remove(obs))
	}

}
