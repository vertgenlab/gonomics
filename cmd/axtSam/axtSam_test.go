package main

import (
	"fmt"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
)

var sequence []dna.Base = dna.StringToBases("GCTTGAAGAGGAGCTGTAGATCAGGCGGAAGATTGGGACCCCCCTTCCCT")
var testSamFile string = "testdata/test.sam"

var axtToSamTests = []struct {
	axtFile   string
	samHeader sam.Header
	output    string
}{
	{"testdata/test.axt", sam.GenerateHeader(fasta.ToChromInfo([]fasta.Fasta{{Name: "chrI", Seq: sequence}}), nil, sam.Unsorted, sam.None), testSamFile},
}

func TestAxtToSamCmd(t *testing.T) {
	for _, v := range axtToSamTests {
		axtToSam(v.axtFile, v.samHeader, v.output)
		testSam, _ := sam.Read(testSamFile)
		var expectedSam sam.Sam = sam.Sam{
			QName: "contig_1",
			Flag:  0,
			RName: "chrI",
			Pos:   1,
			MapQ:  255,
			Cigar: cigar.FromString("38=3D9="),
			RNext: "*",
			PNext: 0,
			TLen:  49,
			Seq:   sequence,
			Qual:  "*",
			Extra: fmt.Sprintf("AS:i:%d\tXS:i:%d\tXE:i:%d", 4021, 1, 47),
		}
		if !sam.Equal(testSam[0], expectedSam) {
			t.Errorf("Error: testing samFile=%v and expectedSam=%v are not equal...", testSam[0], expectedSam)
		}
		os.Remove(testSamFile)
	}
}
