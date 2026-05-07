// Package simulate contains functions for simulation of genomic data with ils
// variants, reads, or regions.
package simulate

import (
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/genePred"
	"github.com/vertgenlab/gonomics/numbers"
	"gonum.org/v1/gonum/mat"
)

// convert a 2D square transition matrix (a, b, c, ...) to (a, a+b, a+b+c, ...) per row
func probRange(transMat *mat.Dense) {
	r, c := transMat.Dims()
	if r != c {
		log.Fatal("Must provide square transition matrix")
	}

	out := mat.NewDense(r, r, nil)
	for i := range r {
		rowSum := 0
		colSum := 0
		for j := range c { // r and c are the same value, named this way for clarity
			rowSum += transMat.At(i, j) // sum at row i
			colSum += transMat.At(j, i) // sum at col i
			out.Set(i, j, rowSum)
			out.Set(j, i, colSum)
		}

		if rowSum != 1{
			log.Fatalf("Must provide transition matrix that sums to 1 across rows, row index %d sums to %d", j, rowSum)
		}
		if colSum != 1{
			log.Fatalf("Must provide transition matrix that sums to 1 across rows, row index %d sums to %d", i, colSum)
		}
	}

	return out
}

// best file practice to read files in cmd
// package deal with data structures, unless for IO
func SimulateIls(roots []*expandedTree.ETree, transMat [][]float64, totalLength int, seed int, outSpecies string, GC float64) []fasta.Fasta{} {
	n := len(roots)
	r, c := transMat.Dims()
	if r != n || c != n {
		log.Fatal("Must provide square transition matrix that matches number of provided phylogenies")
	}
	transMatConverted := probRange(transMat, n)
	
	rand.Seed(seed)

	// initialise sequence of length totalLength
	// random pick start state
	// iterate through and generate next base b ased on transition matrix
	// need to keep track of current state
	
	/// generate an ancestral seq of length totalLength
	/// needs to randomly generate a sequence --> this is only available in command
	// can I reference a function made in command package?
	
	// expandedTree.ReadNewick(filename) (*ETree, error) reades just newick
	anc := fasta.Fasta{Name: "Anc", Seq: RandIntergenicSeq(GC, totalLength)}

	// list of fastas forward simulated from anc
	// if there are S species
	// and N topologies
	// we get S*N total sequences ie a list of length N, each with S sequences
	fastas := make([][]fasta.Fasta, n)
	leafFastas := make([][]fasta.Fasta, n)

	var nodes []*ETree
	for topologyIdx, root := range roots {
		// Simulate takes in a filename, need to save the ancestral sequence??

		func Simulate(randSeqFilename string, root, gene string, deletions bool)
		// now the sequences for each fasta will be in the input tree

		nodes := expandedTree.GetTree(root)
		for i:= 0; i < len(nodes); i++ {
			fastas[topologyIdx] = append(fastas[topologyIdx], *nodes[i].Fasta)
			if nodes[i].Left == nil && nodes[i].Right == nil {
				leafFastas[topologyIdx] = append(leafFastas[topologyIdx], *nodes[i].Fasta)
			}
		}
	}


	// okay, now I have all forward simulations for all the N topologies
	// need to create the final "observed" leaf sequence (equiv human) that takes from each of the N topologies
	// if outSpecies == "", then defaults to the innermost leaf???
	// need to get the "child" sequence name from the tree topo
	speciesIdx := len(fastas[0]) // the innermost species should be the one at the end??? or maybe not, does it 
	// depend on the newick tree order?
	if outSpecies == "" {
		outSpecies = ///// blahblahblah get from the tree
	} else {
		// check if outSpecies is in the tree and get the speciesIdx
	}
	out := fasta.Fasta{Name: outSpecies, Seq: make([]dna.Base, totalLength)}

	currState := rand.Intn(n) // i in matrix
	out.Seq[0] = currState
	transMatConverted := probRange(transitionMat)

	var pickNext float64
	for idx := range totalLength {
		pickNext = rand.Float64() 
		for nextState := range n {
			if pickNext <= transMatConverted.At(currState, nextState) {
				out.Seq[0] = fastas[nextState][speciesIdx].Seq[idx] /// have to also select the correct species
				currState = nextState
			}
		}
	}

	return anc, leaves, out

}
