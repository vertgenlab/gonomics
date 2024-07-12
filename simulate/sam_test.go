package simulate

import (
	"math/rand"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
)

func TestSam(t *testing.T) {
	output := "testdata/actual.sam"

	seed := rand.New(rand.NewSource(1))
	ref := fasta.Read("testdata/eng.fa")
	out := fileio.EasyCreate(output)
	var bw *sam.BamWriter
	header := sam.GenerateHeader(fasta.ToChromInfo(ref), nil, sam.Unsorted, sam.None)
	sam.WriteHeaderToFileHandle(out, header)
	IlluminaPairedSam(ref[0].Name, ref[0].Seq, 100, 150, 500, 50, 0, 0, numbers.BinomialAlias{}, numbers.BinomialAlias{}, 0, out, bw, false, []int{}, seed)

	exception.PanicOnErr(out.Close())
	if !fileio.AreEqual(output, "testdata/expected.sam") {
		t.Errorf("Error in Sam simulation.")
	} else {
		fileio.EasyRemove(output)
	}
}
