package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
	"fmt"
	"log"
)

var seqOneA = dna.StringToBases("ACGTACGTCATCATCATTACTACTAC")
var seqOneB = dna.StringToBases("ACGTACGT")
var seqOneC = dna.StringToBases("ACGTACGTACGTT")
var readWriteTests = []struct {
	filename string // input
	data     []*Node
}{
	{"testdata/testOne.sg", []*Node{{0, seqOneA}, {1, seqOneB}, {2, seqOneC}}},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		actual := Read(test.filename)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
	}
}

func TestWriteAndRead(t *testing.T) {
	var actual []*Node
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		Write(tempFile, test.data)
		actual = Read(tempFile)
		if !AllAreEqual(test.data, actual) {
			t.Errorf("The %s file was not read correctly.", test.filename)
		}
		err := os.Remove(tempFile)
		if err != nil {
			t.Errorf("Deleting temp file %s gave an error.", tempFile)
		}
	}
}

func TestSeeding(t *testing.T) {
        var tileSize int = 25

	genome := Read("testdata/bigGenome.sg")
	reads := fasta.Read("testdata/bigReads.fa")

	log.Printf("Starting to index genome...\n")
        chromPosHash := indexGenome(genome, tileSize)
	log.Printf("Finished indexing!\n")

	log.Printf("Starting to map reads...\n")
        for readNumber := 0; readNumber < len(reads); readNumber++ {
                for tileStart := 0; tileStart < len(reads[readNumber].Seq)-tileSize+1; tileStart++ {
                        codedPositions := chromPosHash[dnaToNumber(reads[readNumber].Seq, tileStart, tileStart+tileSize)]
                        if codedPositions == nil {
                                fmt.Printf("readNumber:%d tileStart:%d notfound\n", readNumber, tileStart)
                        }
                        for hitIndex := 0; hitIndex < len(codedPositions); hitIndex++ {
                                chrom, pos := numberToChromAndPos(codedPositions[hitIndex])
                                fmt.Printf("readNumber:%d tileStart:%d chrom:%d pos:%d\n", readNumber, tileStart, chrom, pos)
                        }
                }
        }
	log.Printf("Finished mapping reads!\n")
}

