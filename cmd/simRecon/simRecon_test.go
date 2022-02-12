package main

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
	"testing"
)

func accEqual(a string, b string) bool {
	var ok = false
	fileA := fileio.Read(a)
	fileB := fileio.Read(b)

	for A := range fileA {
		wordsA := strings.Split(fileA[A], "\t")
		for B := range fileB {
			wordsB := strings.Split(fileB[B], "\t")
			if strings.Compare(wordsA[0], wordsB[0]) == 0 {
				if strings.Compare(wordsA[1], wordsB[1]) == 0 {
					ok = true
				}
			}
		}
	}
	return ok
}

var (
	root     = "testdata/debug.fasta"
	newick   = "testdata/newickShortBranches.txt"
	gene     = "testdata/debug.gp"
	simOut   = "testdata/simOut.fasta"
	leafOut  = "testdata/leafOut.fasta"
	reconOut = "testdata/reconOut.fasta"
	accOut   = "testdata/accOut.txt"
	baseAcc  = "testdata/baseAcc.txt"
)

func TestSimRecon(t *testing.T) {
	SimRecon(root, newick, gene, simOut, leafOut, reconOut, accOut, baseAcc)
	simOk := fasta.AllAreEqualIgnoreOrder(fasta.Read(simOut), fasta.Read("testdata/simOutT.fasta"))
	leafOk := fasta.AllAreEqualIgnoreOrder(fasta.Read(leafOut), fasta.Read("testdata/leafOutT.fasta"))
	reconOk := fasta.AllAreEqualIgnoreOrder(fasta.Read(reconOut), fasta.Read("testdata/reconOutT.fasta"))
	accOk := accEqual(accOut, "testdata/accOutT.txt")
	baseAccOk := accEqual(baseAcc, "testdata/baseAccOutT.txt")
	if simOk && leafOk && reconOk && accOk && baseAccOk != true {
		t.Errorf("one or more testing files did not match expected values, SimMatch: %v, LeafMatch: %v, ReconMatch: %v, AccuracyMatch: %v, BaseAccuracyMatch: %v", simOk, leafOk, reconOk, accOk, baseAccOk)
	}
	fileio.EasyRemove("testdata/simOut.fasta")
	fileio.EasyRemove("testdata/leafOut.fasta")
	fileio.EasyRemove("testdata/reconOut.fasta")
	fileio.EasyRemove("testdata/accOut.txt")
	fileio.EasyRemove("testdata/baseAcc.txt")
}
