package gaf

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestReadAndWrite(t *testing.T) {
	var err error
	actual, header := Read("testdata/test.gaf")
	Write("testdata/out.gaf", actual, header)
	if !fileio.AreEqual("testdata/out.gaf", "testdata/test.gaf") {
		t.Errorf("Error: output is not as expected.")
	} else {
		err = os.Remove("testdata/out.gaf")
		exception.PanicOnErr(err)
	}

	actualChan, headerFromChan := GoReadToChan("testdata/test.gaf")
	out := fileio.EasyCreate("testdata/out.chan.gaf")
	WriteHeaderToFileHandle(out, headerFromChan)
	for i := range actualChan {
		WriteGaf(out, i)
	}
	err = out.Close()
	exception.PanicOnErr(err)
	if !fileio.AreEqual("testdata/out.chan.gaf", "testdata/test.gaf") {
		t.Errorf("Error: output is not as expected with channels.")
	} else {
		err = os.Remove("testdata/out.chan.gaf")
		exception.PanicOnErr(err)
	}
}
