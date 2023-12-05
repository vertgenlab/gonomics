package wig

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ReadWholeGenomeTests = []struct {
	InFile        string
	ChromSizeFile string
	DefaultValue  float64
	OutFile       string
	ExpectedFile  string
}{
	{InFile: "testdata/wholeGenome.wig",
		ChromSizeFile: "testdata/myGenome.chrom.sizes",
		DefaultValue:  0,
		OutFile:       "testdata/test.wholeGenome.wig",
		ExpectedFile:  "testdata/expected.wholeGenome.wig",
	},
}

func TestReadWholeGenome(t *testing.T) {
	var err error
	var wigs map[string]Wig
	for _, v := range ReadWholeGenomeTests {
		wigs = ReadWholeGenome(v.InFile, v.ChromSizeFile, v.DefaultValue)
		WriteMap(v.OutFile, wigs)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error: wig package WriteMap and ReadWholeGenome. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
