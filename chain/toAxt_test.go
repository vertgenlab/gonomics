package chain

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"testing"
)

//TODO: come up with a better test to check conversion. It is possible to be 1 off for keeping track of indices
func TestConvertAxt(t *testing.T) {
	chainfile, _ := Read("testdata/axtTest.chain")
	SortByCoordinates(chainfile, true)
	target := fasta.FastaMap(fasta.Read("testdata/target.fa"))
	query := fasta.FastaMap(fasta.Read("testdata/query.fa"))
	for i := 0; i < len(chainfile); i++ {
		log.Print("chain format:\n")
		log.Printf("%s\n", ChainToString(chainfile[i]))
		log.Print("axt format:\n")
		log.Printf("%s\n", axt.ToString(ChainToAxt(chainfile[i], target[chainfile[i].TName], query[chainfile[i].QName]), i))
	}
}
