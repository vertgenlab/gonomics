// Package simulate contains functions for simulation of genomic data with ils
// variants, reads, or regions.
package simulate

import (
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"gonum.org/v1/gonum/mat"
)

// convert a 2D square transition matrix (a, b, c, ...) to (a, a+b, a+b+c, ...) per row
func probRange(transMat *mat.Dense) *mat.Dense {
	r, c := transMat.Dims()
	if r != c {
		log.Fatal("Must provide square transition matrix")
	}

	out := mat.NewDense(r, r, nil)
	for i := range r {
		rowSum := float64(0)
		colSum := float64(0)
		for j := range c { // r and c are the same value, named this way for clarity
			rowSum += transMat.At(i, j) // sum at row i
			colSum += transMat.At(j, i) // sum at col i
			out.Set(i, j, rowSum)
			out.Set(j, i, colSum)
		}

		if rowSum != 1 {
			log.Fatalf("Must provide transition matrix that sums to 1 across rows, row index %d sums to %d", j, rowSum)
		}
		if colSum != 1 {
			log.Fatalf("Must provide transition matrix that sums to 1 across rows, row index %d sums to %d", i, colSum)
		}
	}

	return out
}

func SimulateIls(roots []*expandedTree.ETree, transMat *mat.Dense, totalLength int, seed int64, outSpecies string, GC float64, gene string, deletions bool) (fasta.Fasta, [][]fasta.Fasta, fasta.Fasta) {
	n := len(roots)
	r, c := transMat.Dims()
	if r != n || c != n {
		log.Fatal("Must provide square transition matrix that matches number of provided phylogenies")
	}
	transMatConverted := probRange(transMat)

	rand.Seed(seed)

	// initialise sequence of length totalLength
	////// random pick start state
	////// iterate through and generate next base b ased on transition matrix
	////// need to keep track of current state

	/// generate an ancestral seq of length totalLength
	anc := []fasta.Fasta{{Name: "Anc", Seq: RandIntergenicSeq(GC, totalLength)}}

	// expandedTree.ReadNewick(filename) (*ETree, error) reades just newick

	// list of fastas forward simulated from anc
	// if there are S species and N topologies,
	//////	we get S*N total sequences ie a list of length N, each with S sequences
	fastas := make([][]fasta.Fasta, n)
	leafFastas := make([][]fasta.Fasta, n)

	var nodes []*expandedTree.ETree

	// need to define gene and deletions
	for topologyIdx, root := range roots {
		SimulateFromSeq(anc, root, gene, deletions)

		nodes = expandedTree.GetTree(root)
		for i := 0; i < len(nodes); i++ {
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
	speciesIdx := -1 // the innermost species should be the one at the end??? or maybe not, does it
	// depend on the newick tree order?
	if outSpecies == "" {
		speciesIdx = len(leafFastas[0]) - 1
		outSpecies = leafFastas[0][speciesIdx].Name
	} else {
		for i, leaf := range leafFastas[0] {
			if leaf.Name == outSpecies {
				speciesIdx = i
				break
			}
		}
	}

	if speciesIdx == -1 {
		log.Fatalf("Could not find species %q in tree", outSpecies)
	}

	ilsBeds := make([][]bed.Bed, n)
	out := fasta.Fasta{Name: outSpecies, Seq: make([]dna.Base, totalLength)}

	prevState := rand.Intn(n) // 0 <= prevState < n
	ilsBeds[prevState] = append(ilsBeds[prevState], bed.Bed{Chrom: "Chr", ChromStart: 0, ChromEnd: 0 + 1, Name: "Sim", FieldsInitialized: 4})
	out.Seq[0] = leafFastas[prevState][speciesIdx].Seq[0] /// need to convert this to a dna.Base

	var pickCurrentStateProb float64
	for idx := range totalLength {
		pickCurrentStateProb = rand.Float64()

		for currentState := range n {
			if pickCurrentStateProb <= transMatConverted.At(prevState, currentState) {
				out.Seq[idx] = leafFastas[currentState][speciesIdx].Seq[idx] /// have to also select the correct species

				if len(ilsBeds[currentState]) == 0 { //// start new bed region run, no existing regions in bed
					ilsBeds[currentState] = append(ilsBeds[currentState], bed.Bed{Chrom: "Chr", ChromStart: idx, ChromEnd: idx + 1, Name: "Sim", FieldsInitialized: 4})
				} else {
					///// what is last region idx
					if ilsBeds[currentState][lastRegionIdx].GetChromEnd() == idx { //// extend existing bed region by 1
						ilsBeds[currentState][lastRegionIdx].ChromEnd = idx + 1
					} else { //// start new bed region run
						ilsBeds[currentState] = append(ilsBeds[currentState], bed.Bed{Chrom: "Chr", ChromStart: idx, ChromEnd: idx + 1, Name: "Sim", FieldsInitialized: 4})
					}
				}

				prevState = currentState

				break
			}
		}

	}

	return anc, leafFastas, out

}
