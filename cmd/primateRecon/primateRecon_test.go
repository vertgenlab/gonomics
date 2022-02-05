package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var PrimateReconTests = []struct {
	Infile       string
	Outfile      string
	ExpectedFile string
	MessyToN     bool
}{
	{"testdata/in.fa", "testdata/tmp.fa", "testdata/expected.fa", true},
}

func TestPrimateRecon(t *testing.T) {
	var err error
	for _, v := range PrimateReconTests {
		primateRecon(v.Infile, v.Outfile, v.MessyToN)
		if !fileio.AreEqual(v.Outfile, v.ExpectedFile) {
			t.Errorf("Error in primateRecon. Output did not match expected.")
		} else {
			err = os.Remove(v.Outfile)
			exception.PanicOnErr(err)
		}
	}
}

// primateReconMle(inFastaFilename string, inTreeFilename string, humanBias bool, probThreshold float64, nonHumanProbThreshold float64, outputFastaFilename string)

func TestAllPossibleOneHuman(t *testing.T) {
	var hum fasta.Fasta = fasta.Fasta{Name: "hg38"}
	var chi fasta.Fasta = fasta.Fasta{Name: "panTro6"}
	var bon fasta.Fasta = fasta.Fasta{Name: "panPan2"}
	var gor fasta.Fasta = fasta.Fasta{Name: "gorGor5"}
	var ora fasta.Fasta = fasta.Fasta{Name: "ponAbe3"}
	var species []fasta.Fasta
	var h, c, b, g, o dna.Base // the values of the human, humanAlt, chimp, ... bases
	h = dna.A                  // human is fixed
	var possibleBases []dna.Base = []dna.Base{dna.A, dna.C, dna.G, dna.T, dna.N, dna.Gap}

	// given that human is A or N or Gap, we will now go through all possible combinations
	for _, c = range possibleBases {
		for _, b = range possibleBases {
			for _, g = range possibleBases {
				for _, o = range possibleBases {
					hum.Seq = append(hum.Seq, h)
					chi.Seq = append(chi.Seq, c)
					bon.Seq = append(bon.Seq, b)
					gor.Seq = append(gor.Seq, g)
					ora.Seq = append(ora.Seq, o)
				}
			}
		}
	}

	species = []fasta.Fasta{hum, chi, bon, gor, ora}
	fasta.Write("testdata/allPossible.fa", species)
}

func TestAllPossibleTwoHumans(t *testing.T) {
	var hum fasta.Fasta = fasta.Fasta{Name: "hg38"}
	var humAlt fasta.Fasta = fasta.Fasta{Name: "hg38alt"}
	var chi fasta.Fasta = fasta.Fasta{Name: "panTro6"}
	var bon fasta.Fasta = fasta.Fasta{Name: "panPan2"}
	var gor fasta.Fasta = fasta.Fasta{Name: "gorGor5"}
	var ora fasta.Fasta = fasta.Fasta{Name: "ponAbe3"}
	var species []fasta.Fasta
	var h, hA, c, b, g, o dna.Base // the values of the human, humanAlt, chimp, ... bases
	h = dna.A                      // human is fixed
	hA = dna.C                     // humanAlt is fixed
	var possibleBases []dna.Base = []dna.Base{dna.A, dna.C, dna.G, dna.T, dna.N, dna.Gap}

	// given that human is A and human alt is C, we will now go through all possible combinations
	for _, c = range possibleBases {
		for _, b = range possibleBases {
			for _, g = range possibleBases {
				for _, o = range possibleBases {
					hum.Seq = append(hum.Seq, h)
					humAlt.Seq = append(hum.Seq, hA)
					chi.Seq = append(chi.Seq, c)
					bon.Seq = append(bon.Seq, b)
					gor.Seq = append(gor.Seq, g)
					ora.Seq = append(ora.Seq, o)
				}
			}
		}
	}

	species = []fasta.Fasta{hum, humAlt, chi, bon, gor, ora}
	fasta.Write("testdata/allPossible.fa", species)
}
