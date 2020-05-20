package chain

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"testing"
) // screen -r 23728.2 41-50

//TODO: come up with a better test to check conversion. It is possible to be 1 off for keeping track of indices
func TestConvertAxt(t *testing.T) {
	chainfile, _ := Read("testdata/axtTest.chain")
	SortByCoordinates(chainfile, true)
	target := fasta.FastaMap(fasta.Read("testdata/target.fa"))
	query := fasta.FastaMap(fasta.Read("testdata/query.fa"))
	for i := 0; i < len(chainfile); i++ {

		log.Printf("%s\n", ChainToString(chainfile[i]))
		log.Printf("%s\n", axt.ToString(ChainToAxt(chainfile[i], target[chainfile[i].TName], query[chainfile[i].QName]), i))
	}
}

func TestConvertChain(t *testing.T) {
	log.SetFlags(0)
	var toyAxt *axt.Axt = &axt.Axt{
		RName:      "chrI",
		RStart:     551,
		REnd:       601,
		QName:      "contig_12",
		QStart:     1,
		QEnd:       51,
		QStrandPos: true,
		Score:      4766,
		RSeq:       dna.StringToBases("CCAGGCTCTATTTTCGGGATTATCCCAAAGGTAGAATGGGTC--GTCACA"),
		QSeq:       dna.StringToBases("CCAGGCTCTATTTTCGGGATTATCCCAAAG--GTAGAATGGGTCGTCACA"),
	}
	var index int = 0
	log.Printf("%s", axt.ToString(toyAxt, 5))
	matches := getChainCounts(toyAxt.RSeq, toyAxt.QSeq)
	index += matches

	log.Printf("Sequence:\n%s\n%s\n", dna.BasesToString(toyAxt.RSeq[:index]), dna.BasesToString(toyAxt.QSeq[:index]))
	target, query := calcMissingBases(toyAxt.RSeq[index:], toyAxt.QSeq[index:])

	log.Printf("Block stats are: %d\t%d\t%d\n", matches, target, query)
	index += common.Max(target, query)
	matches = getChainCounts(toyAxt.RSeq[index:], toyAxt.QSeq[index:])

	log.Printf("the next match is:  length:%d\n%s\n%s\n", matches, dna.BasesToString(toyAxt.RSeq[index:index+matches]), dna.BasesToString(toyAxt.QSeq[index:index+matches]))
	index += matches
	target, query = calcMissingBases(toyAxt.RSeq[index:], toyAxt.QSeq[index:])
	index += common.Max(target, query)
	log.Printf("Block stats are: %d\t%d\t%d\n", matches, target, query)
	log.Printf("curr index = %d, len RSeq=%d", index, len(toyAxt.RSeq))
	log.Printf("Last match is: %d\n\n", matches)

	log.Printf("Put everything together: AxtToChain:\n\n")
	log.Printf("%s", axt.ToString(toyAxt, 5))
	log.Printf("%s\n", ChainToString(AxtToChain(toyAxt, 5)))

}
