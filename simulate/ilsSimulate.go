// Package simulate contains functions for simulation of genomic data with ils
// variants, reads, or regions.
package simulate

import (
	"fmt"
	"log"
	"math"
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
		for j := range c { // r and c are the same value, named this way for clarity
			rowSum += transMat.At(i, j) // sum at row i
			out.Set(i, j, rowSum)
		}

		epsilon := 1e-5
		if math.Abs(rowSum-1.0) > epsilon {
			log.Fatalf("Transition matrix row %d must sum to 1 (within %.1e), got %.10f", i, epsilon, rowSum)
		}
	}

	return out
}

// SimulateILS takes a fasta sequence that is the starting sequence at the root nodes, a pointer to a phylogenetic tree,
// a genePred filename related to the starting sequence, and if deletions should be allowed along with substitutions.
// The starting sequence will then be evolved according to the neutral tree provided and each node in the tree, using
// incomplete lineage separation
// First provided topology should be the non-ILS informed topology
// no gaps
func SimulateIls(roots []*expandedTree.ETree, transMat *mat.Dense, totalLength int, seed int64, chromName string, GC float64, genePred string, deletions bool, leafFastasOnly bool) ([]fasta.Fasta, [][]fasta.Fasta, []bed.Bed, []fasta.Fasta) {
	n := len(roots)
	r, c := transMat.Dims()
	if r != n || c != n {
		log.Fatal("Must provide square transition matrix that matches number of provided phylogenies")
	}
	transMatConverted := probRange(transMat)

	rand.New(rand.NewSource(seed))

	// Randomly generate ancestral sequence of length totalLength
	anc := []fasta.Fasta{{Name: "Anc", Seq: RandIntergenicSeq(GC, totalLength)}}

	// (Repeat for all n topologies) Forward evolve the ancestral sequence.
	// Produces an output fasta with C sequences, one for each node in the topology
	// if there are S species and N topologies,
	// we get S*N total sequences ie a list of length N, each with S sequences
	forwardEvolvedSeqs := make([][]fasta.Fasta, n)

	var nodes []*expandedTree.ETree

	// need to define gene and deletions
	for topologyIdx, root := range roots {
		SimulateFromSeq(anc, root, genePred, deletions)
		nodes = expandedTree.GetTree(root)
		for i := 0; i < len(nodes); i++ {

			if leafFastasOnly {
				if nodes[i].Left == nil && nodes[i].Right == nil {
					forwardEvolvedSeqs[topologyIdx] = append(forwardEvolvedSeqs[topologyIdx], *nodes[i].Fasta)
				}
			} else {
				forwardEvolvedSeqs[topologyIdx] = append(forwardEvolvedSeqs[topologyIdx], *nodes[i].Fasta)
			}
		}
	}

	statePath := GenerateIlsStatePath(totalLength, n, transMatConverted)
	topologyRecord, ilsEvolvedSeqs := CombineIlsSeqs(forwardEvolvedSeqs, statePath, "ilsState")

	return anc, forwardEvolvedSeqs, topologyRecord, ilsEvolvedSeqs

}

func pickState(prevState int, n int, transMatConverted mat.Matrix) int {
	prob := rand.Float64()
	for currentState := range n {
		if prob <= transMatConverted.At(prevState, currentState) {
			return currentState
		}
	}

	return -1
}

func GenerateIlsStatePath(totalLength int, n int, transMatConverted mat.Matrix) []int {
	statePath := make([]int, totalLength)

	prevState := rand.Intn(n)
	statePath[0] = prevState

	for seqIdx := 1; seqIdx < totalLength; seqIdx++ {
		currentState := pickState(prevState, n, transMatConverted)
		if currentState == -1 {
			log.Fatal("Failed to pick next state")
		}

		statePath[seqIdx] = currentState
		prevState = currentState
	}

	return statePath
}

// Generate the topology state record by randomly selecting a starting input state, and then selecting a transition
// Combine the n forward-evolved multifastas using the transition matrix
// Have to iterate over all C sequences in the multifasta!!!!!
func CombineIlsSeqs(forwardEvolvedSeqs [][]fasta.Fasta, states []int, chromName string) ([]bed.Bed, []fasta.Fasta) {
	numTopologies := len(forwardEvolvedSeqs)
	numSpecies := len(forwardEvolvedSeqs[0])
	totalLength := len(forwardEvolvedSeqs[0][0].Seq)

	if len(states) != totalLength {
		log.Fatalf("state path length %d does not match sequence length %d", len(states), totalLength)
	}

	speciesIndex := make([]map[string]int, numTopologies)
	for topoIdx := range speciesIndex {
		speciesIndex[topoIdx] = make(map[string]int)
		for seqIdx := range forwardEvolvedSeqs[topoIdx] {
			speciesIndex[topoIdx][forwardEvolvedSeqs[topoIdx][seqIdx].Name] = seqIdx
		}
	}

	recordName := make(map[int]string)
	for i := 0; i < numTopologies; i++ {
		recordName[i] = fmt.Sprintf("V%d", i)
	}

	ilsEvolvedSeqs := make([]fasta.Fasta, numSpecies)
	for speciesIdx := 0; speciesIdx < numSpecies; speciesIdx++ {
		ilsEvolvedSeqs[speciesIdx].Name = forwardEvolvedSeqs[0][speciesIdx].Name
		ilsEvolvedSeqs[speciesIdx].Seq = make([]dna.Base, totalLength)
	}

	var topologyRecord []bed.Bed
	var currentState int
	var currentName string
	var srcIdx int

	for seqIdx := 0; seqIdx < totalLength; seqIdx++ {
		currentState = states[seqIdx]
		currentName = recordName[currentState]

		if len(topologyRecord) > 0 && topologyRecord[len(topologyRecord)-1].ChromEnd == seqIdx && topologyRecord[len(topologyRecord)-1].Name == currentName {
			topologyRecord[len(topologyRecord)-1].ChromEnd = seqIdx + 1
		} else {
			topologyRecord = append(topologyRecord, bed.Bed{Chrom: chromName, ChromStart: seqIdx, ChromEnd: seqIdx + 1, Name: currentName, FieldsInitialized: 4})
		}

		for species, specIdx := range speciesIndex[0] {
			srcIdx = speciesIndex[currentState][species]
			ilsEvolvedSeqs[specIdx].Seq[seqIdx] = forwardEvolvedSeqs[currentState][srcIdx].Seq[seqIdx]
		}
	}

	return topologyRecord, ilsEvolvedSeqs
}
