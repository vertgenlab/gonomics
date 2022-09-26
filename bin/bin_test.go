package bin

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"testing"
)

var (
	genome1     = "testdata/testBreakBinning.fa"
	expBases1   = dna.StringToBases("ACGTA")
	expBases2   = dna.StringToBases("CGT")
	expBases3   = dna.StringToBases("AC")
	expBases4   = dna.StringToBases("GTACG")
	expBases5   = dna.StringToBases("ACGTA")
	expBases6   = dna.StringToBases("AC")
	expBases7   = dna.StringToBases("AC")
	expNames1   = "chr1"
	expNames2   = "chr2"
	expNames3   = "chr3"
	expNames4   = "chr4"
	expNames5   = "chr5"
	genome2     = "testdata/testContigs.fa"
	expNoBreak1 = fasta.Fasta{"chr1", dna.StringToBases("ACGTACGT")}
	expNoBreak2 = fasta.Fasta{"chr2", dna.StringToBases("ACGT")}
	expNoBreak3 = fasta.Fasta{"chr3", dna.StringToBases("ACG")}
	expNoBreak4 = fasta.Fasta{"chr4", dna.StringToBases("T")}
)

func TestBinFasta(t *testing.T) {
	gen := fasta.Read(genome1)
	out := BinFasta(gen, 5)

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
	}
}

//expectation: 2 bins, min 4, chr1 ACGTACGT, chr2-4, ACGTACGT, name convention will be chr2_0chr3_4chr4_7.fa
func TestBinGenomeNoBreaks(t *testing.T) {
	records := fasta.Read(genome2)
	bins := BinGenomeNoBreaks(records, 2, -1)
	binsMin := BinGenomeNoBreaks(records, 0, 4)

	log.Print(bins)
	log.Print(binsMin)
	if len(bins) != 2 || len(binsMin) != 2 {
		log.Fatalf("wrong number of bins created. bins: %v, binsMin: %v.", len(bins), len(binsMin))
	}

	for i := range bins {
		value, _ := bins[i]
		minValue, _ := binsMin[i]
		if i == 0 {
			if !fasta.IsEqual(value[0], expNoBreak1) || !fasta.IsEqual(minValue[0], expNoBreak1) {
				log.Fatalf("First fasta in bin: %s %s or first fasta in binsMin: %s %s didn't match expected value: %s %s.", value[0].Name, value[0].Seq, minValue[0].Name, minValue[0].Seq, expNoBreak1.Name, expNoBreak1.Seq)
			}
		} else {
			if !fasta.IsEqual(value[0], expNoBreak2) || !fasta.IsEqual(minValue[0], expNoBreak2) {
				log.Fatalf("Fasta in bin: %s %s or fasta in binsMin: %s %s didn't match expected value: %s %s.", value[0].Name, value[0].Seq, minValue[0].Name, minValue[0].Seq, expNoBreak2.Name, expNoBreak2.Seq)
			} else if !fasta.IsEqual(value[1], expNoBreak3) || !fasta.IsEqual(minValue[1], expNoBreak3) {
				log.Fatalf("Fasta in bin: %s %s or fasta in binsMin: %s %s didn't match expected value: %s %s.", value[1].Name, value[1].Seq, minValue[1].Name, minValue[1].Seq, expNoBreak3.Name, expNoBreak3.Seq)
			} else if !fasta.IsEqual(value[2], expNoBreak4) || !fasta.IsEqual(minValue[2], expNoBreak4) {
				log.Fatalf("Fasta in bin: %s %s or fasta in binsMin: %s %s didn't match expected value: %s %s.", value[2].Name, value[2].Seq, minValue[2].Name, minValue[2].Seq, expNoBreak4.Name, expNoBreak4.Seq)
			}
		}
	}
}
