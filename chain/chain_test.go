package chain

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var dir []string = []string{"testdata/big.chain", "testdata/twoChainZ.chain"}

func TestReadAndWrite(t *testing.T) {
	var tmpFile string
	var err error
	for _, v := range dir {
		tmpFile = v + ".tmp"
		records, header := Read(v)
		Write(tmpFile, records, header)
		if !fileio.AreEqual(tmpFile, v) {
			t.Errorf("Error in chain package read and write.")
		} else {
			err = os.Remove(tmpFile)
			exception.PanicOnErr(err)
		}
	}
}
