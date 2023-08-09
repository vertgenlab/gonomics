package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func TestAllPossible(t *testing.T) {
	var err error
	var mouse fasta.Fasta = fasta.Fasta{Name: "mm10"}
	var rat fasta.Fasta = fasta.Fasta{Name: "rn7"}
	var hamster fasta.Fasta = fasta.Fasta{Name: "criGriChoV2"}
	var squirrel fasta.Fasta = fasta.Fasta{Name: "speTri2"}
	var species []fasta.Fasta
	var m, r, h, s dna.Base //the values of the current base for each species
	m = dna.A               //mouse is fixed
	var possibleBases []dna.Base = []dna.Base{dna.A, dna.C, dna.G, dna.T, dna.N, dna.Gap}

	// given that mouse is A, we will now go through all possible combinations
	for _, r = range possibleBases {
		for _, h = range possibleBases {
			for _, s = range possibleBases {
				mouse.Seq = append(mouse.Seq, m)
				rat.Seq = append(rat.Seq, r)
				hamster.Seq = append(hamster.Seq, h)
				squirrel.Seq = append(squirrel.Seq, s)
			}
		}
	}

	species = []fasta.Fasta{mouse, rat, hamster, squirrel}
	fasta.Write("testdata/allPossible.fa", species)

	mouseReconMraMle("testdata/allPossible.fa", "testdata/test.mraMleMouseBias.fa", "testdata/4d.mod", 0.0, 0.8)
	if !fileio.AreEqual("testdata/test.mraMleMouseBias.fa", "testdata/expected.mraMleMouseBias.fa") {
		t.Errorf("Error in mouseRecon, mouse-biased mleMra reconstruction. Output was not as expected.")
	} else {
		err = os.Remove("testdata/test.mraMleMouseBias.fa")
		exception.PanicOnErr(err)
	}

	mouseReconMraMle("testdata/allPossible.fa", "testdata/test.mraMleMouseBias0.fa", "testdata/4d.mod", 0.0, 0.0)
	if !fileio.AreEqual("testdata/test.mraMleMouseBias0.fa", "testdata/expected.mraMleMouseBias0.fa") {
		t.Errorf("Error in mouseRecon, mouse-biased mleMra reconstruction. Output was not as expected.")
	} else {
		err = os.Remove("testdata/test.mraMleMouseBias0.fa")
		exception.PanicOnErr(err)
		err = os.Remove("testdata/allPossible.fa") //remove allPossible once all tests have passed.
		exception.PanicOnErr(err)
	}
}
