package axt

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"testing"
)

func TestSamFileConvert(t *testing.T) {

	//chr19 3001012 3001075 chr11 70568380 70568443 - 3500
	var testAxt = &Axt{
		RName:      "chr19",
		RStart:     1,
		REnd:       30,
		QName:      "chr11",
		QStart:     2,
		QEnd:       31,
		QStrandPos: false,
		Score:      3500,
		RSeq:       dna.StringToBases("TCAGCTCATAAATCACCTCC----ACAAGC"),
		QSeq:       dna.StringToBases("TCTG--CATAAACCACCTGCCATGACAAGC"),
	}
	samFromAxt := AxtToSam(testAxt)

	var answerSam *sam.SamAln = &sam.SamAln{
		QName: "chr11",
		Flag:  16,
		RName: "chr19",
		Pos:   1,
		MapQ:  255, // mapping quality setting to 255 because we are not calculating it
		Cigar: cigar.FromString("2=1X1=2D6=1X5=1X1=4I6="),
		RNext: "*",
		PNext: 0,
		TLen:  29, //Could leave at zero or make TLen be the length of alignment, start and end (not sure if i can get target length from an axt)
		Seq:   dna.StringToBases("TCTGCATAAACCACCTGCCATGACAAGC"),
		Qual:  "*",
		Extra: fmt.Sprintf("AS:i:%d\tXS:i:%d\tXE:i:%d", 3500, 2, 31),
	}
	if !sam.IsEqual(samFromAxt, answerSam) {
		t.Errorf("Error: Axt to sam is not converting the correct sam file...\n")
	}

}
