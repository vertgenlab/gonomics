package starrSeq

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

func TestPlasmidUMI(t *testing.T) {
	inBed := "testdata/in.umi.bed"
	inSam := "testdata/in.umi.sam"
	out := "testdata/out.umi.txt"
	exp := "testdata/exp.umi.txt"
	PlasmidUMI(inSam, inBed, out)
	if !fileio.AreEqual(exp, out) {
		t.Errorf("expected and observed files are not equal. exp: %s, obs: %s", exp, out)
	} else {
		err := os.Remove(out)
		exception.PanicOnErr(err)
	}
}
