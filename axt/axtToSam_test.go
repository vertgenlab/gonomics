package axt

import (

	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"testing"
)

func TestSamFileConvert(t *testing.T) {
	//chr19 3001012 3001075 chr11 70568380 70568443 - 3500
	var testAxt  = &Axt{
		RName:      "chr19",
		RStart:     3001012,
		REnd:       3001075,
		QName:      "chr11",
		QStart:     70568380,
		QEnd:       70568443,
		QStrandPos: false,
		Score:      3500,
		RSeq:       dna.StringToBases("TCAGCTCATAAATCACCTCCTGCCACAAGCCTGGCCTGGTCCCAGGAGAGTGTCCAGGCTCAGA"),
		QSeq:       dna.StringToBases("TCTGTTCATAAACCACCTGCCATGACAAGCCTGGCCTGTTCCCAAGACAATGTCCAGGCTCAGA"),
	}
	log.Printf("%s\n", sam.SamAlnToString(AxtToSam(testAxt)))
	log.Printf("Toy example: PASS!")
	testFile := fileio.EasyOpen("testdata/rabsToGasAcu1_filteredByScore.axt")
	reader := make(chan *Axt)
	defer testFile.Close()
	go ReadToChan(testFile, reader)
	for each := range reader {
		log.Printf("%s\n", sam.SamAlnToString(AxtToSam(each)))
	}

}

