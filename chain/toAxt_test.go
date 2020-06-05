package chain

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"testing"
)

func TestConvertAxt(t *testing.T) {
	chainfile, _ := Read("testdata/axtTest.chain")
	SortByCoordinates(chainfile, true)
	target := fasta.FastaMap(fasta.Read("testdata/target.fa"))
	query := fasta.FastaMap(fasta.Read("testdata/query.fa"))
	for i := 0; i < len(chainfile); i++ {
		log.Print("chain format:\n")
		log.Printf("%s\n", ToString(chainfile[i]))
		log.Print("axt format:\n")
		log.Printf("%s\n", axt.ToString(ChainToAxt(chainfile[i], target[chainfile[i].TName], query[chainfile[i].QName]), i))
	}
}
