package starrSeq

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestMkref(t *testing.T) {
	Mkref("testdata/upstream.fa", "testdata/constructs.fa", "testdata/downstream.fa", "testdata/testOut")
	if !fileio.AreEqual("testdata/testOut.fa", "testdata/exp.mkref.fa") {
		t.Errorf("Error in mkref, files aren't equal: %s, %s", "testOut.fa", "exp.mkref.fa")
	} else {
		err := os.Remove("testdata/testOut.fa")
		exception.PanicOnErr(err)
	}
	if !fileio.AreEqual("testdata/testOut.bed", "testdata/exp.mkref.bed") {
		t.Errorf("Error in mkref, files aren't equal: %s, %s", "testOut.bed", "exp.mkref.bed")
	} else {
		err := os.Remove("testdata/testOut.bed")
		exception.PanicOnErr(err)
	}
}
