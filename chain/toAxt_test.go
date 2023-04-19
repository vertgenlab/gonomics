package chain

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestToAxt(t *testing.T) {
	var err error
	chainfile, _ := Read("testdata/axtTest.chain")
	SortByCoordinates(chainfile, true)
	target := fasta.ToMap(fasta.Read("testdata/target.fa"))
	query := fasta.ToMap(fasta.Read("testdata/query.fa"))
	answer := AllToAxt(chainfile, target, query)
	axt.Write("testdata/tmp.AxtTest.axt", answer)
	if !fileio.AreEqual("testdata/tmp.AxtTest.axt", "testdata/expected.ToAxt.axt") {
		t.Errorf("Error in ToAxt. Output was not as expected.")
	} else {
		err = os.Remove("testdata/tmp.AxtTest.axt")
		exception.PanicOnErr(err)
	}
}
