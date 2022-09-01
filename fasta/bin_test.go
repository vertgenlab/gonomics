package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

var (
	genome    = "testdata/testContigs.fa"
	expBases1 = dna.StringToBases("ACGTA")
	expBases2 = dna.StringToBases("CGT")
	expBases3 = dna.StringToBases("AC")
	expBases4 = dna.StringToBases("GTACG")
	expBases5 = dna.StringToBases("ACGTA")
	expBases6 = dna.StringToBases("AC")
	expBases7 = dna.StringToBases("AC")
	expNames1 = "chr1"
	expNames2 = "chr2"
	expNames3 = "chr3"
	expNames4 = "chr4"
	expNames5 = "chr5"
	outName   = "testdata/bin.fa"
	outNames  = []string{"testdata/bin.0.fa", "testdata/bin.1.fa", "testdata/bin.2.fa", "testdata/bin.3.fa", "testdata/bin.4.fa"}
)

func TestBinFasta(t *testing.T) {
	gen := Read(genome)
	out := BinFasta(gen, 5)

	log.Print("HELLO")

	for m := range out {
		if m == 0 {
			if out[m][0].Name != expNames1 {
				t.Fatalf("expected first map entry to have name %s but showed %s instead", expNames1, out[m][0].Name)
			} else if dna.CompareSeqsCaseSensitive(expBases1, out[m][0].Seq) != 0 {
				t.Fatalf("expected first map entry sequence to be %s but was %s instead", dna.BasesToString(expBases1), dna.BasesToString(out[m][0].Seq))
			}
		} else if m == 1 {
			if out[m][0].Name != expNames1 || out[m][1].Name != expNames2 {
				t.Fatalf("expected second map entry to have names %s and %s but showed %s and %s instead", expNames1, expNames2, out[m][0].Name, out[m][1].Name)
			} else if dna.CompareSeqsCaseSensitive(expBases2, out[m][0].Seq) != 0 || dna.CompareSeqsCaseSensitive(expBases3, out[m][1].Seq) != 0 {
				t.Fatalf("expected second map entry sequence to be %s but was %s instead", dna.BasesToString(expBases1), dna.BasesToString(out[m][0].Seq))
			}
		} else if m == 2 {
			if out[m][0].Name != expNames2 {
				t.Fatalf("expected third map entry to have name %s but showed %s instead", expNames3, out[m][0].Name)
			} else if dna.CompareSeqsCaseSensitive(expBases4, out[m][0].Seq) != 0 {
				t.Fatalf("expected third map entry sequence to be %s but was %s instead", dna.BasesToString(expBases4), dna.BasesToString(out[m][0].Seq))
			}
		} else if m == 3 {
			if out[m][0].Name != expNames3 {
				t.Fatalf("expected fourth map entry to have name %s but showed %s instead", expNames1, out[m][0].Name)
			} else if dna.CompareSeqsCaseSensitive(expBases5, out[m][0].Seq) != 0 {
				t.Fatalf("expected fourth map entry sequence to be %s but was %s instead", dna.BasesToString(expBases1), dna.BasesToString(out[m][0].Seq))
			}
		} else if m == 4 {
			if out[m][0].Name != expNames4 || out[m][1].Name != expNames5 {
				t.Fatalf("expected fifth map entry to have name %s and %s but showed %s and %s instead", expNames4, expNames5, out[m][0].Name, out[m][1].Name)
			} else if dna.CompareSeqsCaseSensitive(expBases6, out[m][0].Seq) != 0 || dna.CompareSeqsCaseSensitive(expBases7, out[m][1].Seq) != 0 {
				t.Fatalf("expected fifth map entry sequence to be %s and %s but was %s and %s instead", dna.BasesToString(expBases6), dna.BasesToString(expBases7), dna.BasesToString(out[m][0].Seq), dna.BasesToString(out[m][1].Seq))
			}
		} else {
			t.Fatal("Map does not contain expected values")
		}
		dir, file := filepath.Split(outName)
		outFile := strings.TrimSuffix(file, ".fa")
		outFileName := dir + outFile + "." + fileio.IntToString(m) + ".fa"
		log.Print(outFileName)
		Write(outFileName, out[m])
	}
	for r := range outNames {
		os.Remove(outNames[r])
	}
}
