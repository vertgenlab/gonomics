package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

func TestAxtToFasta(t *testing.T) {
	log.SetFlags(0)
	var testAxt = axt.Axt{
		RName:      "chr19",
		RStart:     11,
		REnd:       40,
		QName:      "chr11",
		QStart:     2,
		QEnd:       31,
		QStrandPos: false,
		Score:      3500,
		RSeq:       dna.StringToBases("TCTGNNNNTAAACCACCNNNNATGACAAGC"),
		QSeq:       dna.StringToBases("TCAGCTCATAAATCACCTCCCATGACAAGC"),
	}
	var testFa = &fasta.Fasta{
		Name: "chr19",
		Seq:  dna.StringToBases("AAAAAAAAAATCTGNNNNTAAACCACCNNNNATGACAAGC"),
	}

	seq := axtSeq(testAxt, testFa.Seq)
	log.Printf("%s %s\n", seq.Name, dna.BasesToString(seq.Seq))
}
