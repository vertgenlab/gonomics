package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"log"
	"testing"
	//	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
)

func TestAxtToFasta(t *testing.T) {
	log.SetFlags(0)
	var testAxt = &axt.Axt{
		RName:      "chr19",
		RStart:     1,
		REnd:       30,
		QName:      "chr11",
		QStart:     2,
		QEnd:       31,
		QStrandPos: false,
		Score:      3500,
		RSeq:       dna.StringToBases("TCAGCTCATAAATCACCTCCCATGACAAGC"),
		QSeq:       dna.StringToBases("TCTGNNNNTAAACCACCNNNNATGACAAGC"),
	}
}
