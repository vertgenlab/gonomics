package simulate

import (
	"math/rand"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
)

func TestSam(t *testing.T) {
	var err error
	rand.Seed(1)
	ref := fasta.Read("testdata/eng.fa")
	out := fileio.EasyCreate("testdata/actual.sam")
	var bw *sam.BamWriter
	header := sam.GenerateHeader(fasta.ToChromInfo(ref), nil, sam.Unsorted, sam.None)
	sam.WriteHeaderToFileHandle(out, header)
	IlluminaPairedSam(ref[0].Name, ref[0].Seq, 100, 150, 50, 50, 0, numbers.BinomialAlias{}, out, bw, false)
	err = out.Close()
	exception.PanicOnErr(err)
	if !fileio.AreEqual("testdata/actual.sam", "testdata/expected.sam") {
		t.Errorf("Error in Sam simulation.")
	} else {
		err = os.Remove("testdata/actual.sam")
		exception.PanicOnErr(err)
	}
}
