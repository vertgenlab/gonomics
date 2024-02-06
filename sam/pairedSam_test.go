package sam

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestReadSamPeToChan(t *testing.T) {
	out := "testdata/out.pe.sam"
	exp := "testdata/exp.pe.sam"
	ch, head := GoReadSamPeToChan("testdata/pe.sam")
	o := fileio.EasyCreate(out)
	WriteHeaderToFileHandle(o, head)
	for i := range ch {
		WriteSamPeToFileHandle(o, i)
	}
	err := o.Close()
	exception.PanicOnErr(err)

	if !fileio.AreEqual(exp, out) {
		t.Errorf("Expected and output files are not identical: %s, %s", exp, out)
	} else {
		err = os.Remove(out)
		exception.PanicOnErr(err)
	}
}
