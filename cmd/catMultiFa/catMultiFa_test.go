package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var CatMultiFaTests = []struct {
	FileList   []string
	O          string
	LineLength int
	Expected   string
}{
	{[]string{"testdata/file1.fa", "testdata/file2.fa", "testdata/file3.fa"}, "testdata/tmp.fa", 50, "testdata/expected.fa"},
}

func TestCatMultiFa(t *testing.T) {
	var err error
	for _, v := range CatMultiFaTests {
		catMultiFa(v.FileList, v.O, v.LineLength)
		if !fileio.AreEqual(v.O, v.Expected) {
			t.Errorf("Error in catMultiFa.")
		} else {
			err = os.Remove(v.O)
			exception.PanicOnErr(err)
		}
	}
}
