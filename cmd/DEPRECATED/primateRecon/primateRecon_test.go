package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestAllPossibleOneHumanNonGeneric(t *testing.T) {
	var err error
	var hum fasta.Fasta = fasta.Fasta{Name: "hg38"}
	var chi fasta.Fasta = fasta.Fasta{Name: "panTro6"}
	var bon fasta.Fasta = fasta.Fasta{Name: "panPan2"}
	var gor fasta.Fasta = fasta.Fasta{Name: "gorGor5"}
	var ora fasta.Fasta = fasta.Fasta{Name: "ponAbe3"}
	var species []fasta.Fasta
	var h, c, b, g, o dna.Base // the values of the human, humanAlt, chimp, ... bases
	var hBases []dna.Base = []dna.Base{dna.A, dna.N, dna.Gap}
	var possibleBases []dna.Base = []dna.Base{dna.A, dna.C, dna.G, dna.T, dna.N, dna.Gap}

	// given that human is A or N or Gap, we will now go through all possible combinations
	for _, h = range hBases {
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
	}

	species = []fasta.Fasta{hum, chi, bon, gor, ora}
	fasta.Write("testdata/allPossible.oneHuman.fa", species)

	primateRecon("testdata/allPossible.oneHuman.fa", "testdata/out.humanBiasedParsimony.fa", false)
	if !fileio.AreEqual("testdata/out.humanBiasedParsimony.fa", "testdata/expected.humanBiasedParsimony.fa") {
		t.Errorf("Error in primateRecon, human biased manual algorithm. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.humanBiasedParsimony.fa")
		exception.PanicOnErr(err)
	}

	primateRecon("testdata/allPossible.oneHuman.fa", "testdata/out.ParsimonyMessyToN.fa", true)
	if !fileio.AreEqual("testdata/out.ParsimonyMessyToN.fa", "testdata/expected.ParsimonyMessyToN.fa") {
		t.Errorf("Error in primateRecon, human biased manual algorithm. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.ParsimonyMessyToN.fa")
		exception.PanicOnErr(err)
	}

	primateReconHcaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", true, false, 0.0, 0, false, "testdata/out.humanBiasedMleNoThreshold.fa")
	if !fileio.AreEqual("testdata/out.humanBiasedMleNoThreshold.fa", "testdata/expected.humanBiasedMleNoThreshold.fa") {
		t.Errorf("Error in primateRecon, human biased version nonHumanProb no threshold. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.humanBiasedMleNoThreshold.fa")
		exception.PanicOnErr(err)
	}

	primateReconHcaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", true, false, 0.0, 0.99, false, "testdata/out.humanBiasedMle99.fa")
	if !fileio.AreEqual("testdata/out.humanBiasedMle99.fa", "testdata/expected.humanBiasedMle99.fa") {
		t.Errorf("Error in primateRecon, human biased version nonHumanProb 99. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.humanBiasedMle99.fa")
		exception.PanicOnErr(err)
	}

	primateReconHcaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", true, false, 0.0, 0.8, false, "testdata/out.humanBiasedMle80.fa")
	if !fileio.AreEqual("testdata/out.humanBiasedMle80.fa", "testdata/expected.humanBiasedMle80.fa") {
		t.Errorf("Error in primateRecon, human biased version nonHumanProb 80. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.humanBiasedMle80.fa")
		exception.PanicOnErr(err)
	}
	primateReconHcaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", false, true, 0.0, 0.8, false, "testdata/out.chimpBiasedMle80.fa")
	if !fileio.AreEqual("testdata/out.chimpBiasedMle80.fa", "testdata/expected.chimpBiasedMle80.fa") {
		t.Errorf("Error in primateRecon, chimp biased version nonChimpProb 80. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.chimpBiasedMle80.fa")
		exception.PanicOnErr(err)
	}
	primateReconHgaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", 0.0, 0.8, false, "testdata/out.gorillaBiasedHgaMle80.fa")
	if !fileio.AreEqual("testdata/out.gorillaBiasedHgaMle80.fa", "testdata/expected.gorillaBiasedHgaMle80.fa") {
		t.Errorf("Error in primateRecon, gorilla biased hga version nonBiasProb 80. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.gorillaBiasedHgaMle80.fa")
		exception.PanicOnErr(err)
	}
	primateReconHgaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.mod", 0.0, 0.99, false, "testdata/out.gorillaBiasedHgaMle99.fa")
	if !fileio.AreEqual("testdata/out.gorillaBiasedHgaMle99.fa", "testdata/expected.gorillaBiasedHgaMle99.fa") {
		t.Errorf("Error in primateRecon, gorilla biased hga version nonBiasProb 99. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.gorillaBiasedHgaMle99.fa")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/allPossible.oneHuman.fa") //we delete allPossible only if all the tests have passed.
		exception.PanicOnErr(err)
	}
}

func TestAllPossibleOneHumanGenericNames(t *testing.T) {
	var err error
	var hum fasta.Fasta = fasta.Fasta{Name: "human"}
	var chi fasta.Fasta = fasta.Fasta{Name: "chimp"}
	var bon fasta.Fasta = fasta.Fasta{Name: "bonobo"}
	var gor fasta.Fasta = fasta.Fasta{Name: "gorilla"}
	var ora fasta.Fasta = fasta.Fasta{Name: "orangutan"}
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
	fasta.Write("testdata/allPossible.oneHuman.fa", species)

	primateReconHcaMle("testdata/allPossible.oneHuman.fa", "testdata/4d.genericNames.mod", true, false, 0.0, 0.99, true, "testdata/out.humanBiasedMle99.genericNames.fa")
	if !fileio.AreEqual("testdata/out.humanBiasedMle99.genericNames.fa", "testdata/expected.humanBiasedMle99.genericNames.fa") {
		t.Errorf("Error in primateRecon, human biased version nonHumanProb 99. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.humanBiasedMle99.genericNames.fa")
		exception.PanicOnErr(err)
	}
}

func TestAllPossibleTwoHumans(t *testing.T) {
	var err error
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
					humAlt.Seq = append(humAlt.Seq, hA)
					chi.Seq = append(chi.Seq, c)
					bon.Seq = append(bon.Seq, b)
					gor.Seq = append(gor.Seq, g)
					ora.Seq = append(ora.Seq, o)
				}
			}
		}
	}

	species = []fasta.Fasta{hum, humAlt, chi, bon, gor, ora}
	fasta.Write("testdata/allPossible.twoHumans.fa", species)
	primateReconHcaMle("testdata/allPossible.twoHumans.fa", "testdata/4d.2h.mod", false, false, 0.90, 0.0, false, "testdata/out.unbiased90.fa")
	if !fileio.AreEqual("testdata/out.unbiased90.fa", "testdata/expected.unbiased90.fa") {
		t.Errorf("Error in primateRecon, unbiased 90. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.unbiased90.fa")
		exception.PanicOnErr(err)
	}
	primateReconHcaMle("testdata/allPossible.twoHumans.fa", "testdata/4d.2h.mod", false, false, 0.99, 0.0, false, "testdata/out.unbiased99.fa")
	if !fileio.AreEqual("testdata/out.unbiased99.fa", "testdata/expected.unbiased99.fa") {
		t.Errorf("Error in primateRecon, unbiased 99. Output was not as expected.")
	} else {
		err = os.Remove("testdata/out.unbiased99.fa")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/allPossible.twoHumans.fa")
		exception.PanicOnErr(err)
	}
}
