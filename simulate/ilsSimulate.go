// Package simulate contains functions for simulation of genomic data with ils
// variants, reads, or regions.
package simulate

import (
	"fmt"
	"log"
	"math"
	"math/rand"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers/matrix"
	"github.com/vertgenlab/gonomics/reconstruct"
	"github.com/vertgenlab/gonomics/wig"
	"gonum.org/v1/gonum/mat"
)

// cumulativeTransitionMatrix converts a 2D square transition matrix to a cumulative transition matrix
// for example, given a stochastic matrix where row 1 = (a, b, c, ..., n) such that row 1 sums to 1,
// and a, b, c... describe the probability of the state 1 converting to any state N,
// We then convert the row 1 to (a, a+b, a+b+c, ..., a+b+⋯+n) such that the last value in the row is 1.
func cumulativeTransitionMatrix(transMat *mat.Dense) *mat.Dense {
	r, c := transMat.Dims()
	if r != c {
		log.Fatal("Must provide square transition matrix")
	}

	var curProb float64
	out := mat.NewDense(r, r, nil)
	for i := range r {
		rowSum := float64(0)
		for j := range c { // r and c are the same value, named this way for clarity
			curProb = transMat.At(i, j)
			if curProb < 0 || curProb > 1 {
				log.Fatalf("Cannot provide negative transition probability or probability over 1.")
			}
			rowSum += curProb // sum at row i
			out.Set(i, j, rowSum)
		}

		epsilon := 1e-5
		if math.Abs(rowSum-1.0) > epsilon {
			log.Fatalf("Transition matrix row %d must sum to 1 (within %.1e), got %.10f", i, epsilon, rowSum)
		}
	}

	return out
}

// SimulateIls simulates sequence evolution across a set of phylogenetic trees,
// incorporating incomplete lineage sorting (ILS).
//   - roots: Pointer to a set of phylogenetic trees
//   - transMat: Transition matrix defining transitions between phylogenetic
//     topologies and their unit branch lengths
//   - ancSeq: FASTA sequence used as the ancestral sequence at the root nodes
//     May be empty if totalLength is nonzero and ancName is nonempty
//   - totalLength: Length of the ancestral sequence to generate if ancSeq is empty
//   - ancName: Name of the ancestral sequence; required if ancSeq is empty
//   - seed: Random seed used for sequence evolution
//   - chromName: Chromosome name used to construct output sequence names
//   - leafFastasOnly: If true, outputs only leaf-node sequences
//   - genePred: GenePred filename associated with the starting sequence
//   - deletions: Whether to allow deletions in addition to substitutions
//
// Uses defaultSubstitutionMatrix if no substitution matrix file is provided.
// The ancestral sequence is evolved along the supplied phylogenetic trees,
// accounting for topology changes due to ILS. No gaps are introduced.
//
// Outputs:
//   - ancestral sequence
//   - multi-FASTA file containing forward-evolved sequences without ILS
//   - BED file identifying the topology associated with each position
//   - ILS-informed sequences multifa
func SimulateIls(roots []*expandedTree.ETree, transMat *mat.Dense, ancSeq []fasta.Fasta, totalLength int, seed int64, chromName string, leafFastasOnly bool, substitutionMatrixFile string, unitBranchLength float64, ancName string) ([]fasta.Fasta, [][]fasta.Fasta, []bed.Bed, []fasta.Fasta) {

	// set entire genome to be 1 big gene
	numRoots := len(roots)
	r, c := transMat.Dims()
	if r != numRoots || c != numRoots {
		log.Fatal("Must provide square transition matrix that matches number of provided phylogenies")
	}

	transMatConverted := cumulativeTransitionMatrix(transMat)

	// TODO convert to seed object
	rand.Seed(seed)
	var anc []fasta.Fasta
	if len(ancSeq) == 0 {
		if ancName == "" {
			log.Fatal("Must provide name of ancestral sequence to be generated, if no input ancestral sequence provide.")
		}
		anc = []fasta.Fasta{{Name: ancName, Seq: RandIntergenicSeq(GC, totalLength)}}
	} else {
		// Randomly generate ancestral sequence of length totalLength
		anc = ancSeq
		totalLength = len(anc[0].Seq)
	}

	// (Repeat for all n topologies) Forward evolve the ancestral sequence.
	// Produces an output fasta with C sequences
	// 
	, one for each node in the topology
	// if there are S species and N topologies,
	// we get S*N total sequences ie a list of length N, each with S sequences
	forwardEvolvedSeqs := make([][]fasta.Fasta, numRoots)

	var nodes []*expandedTree.ETree

	// TODO: there should be a check that the name of the ancestral sequence matches the name of the topmost root
	for topologyIdx, root := range roots {
		root.Fasta = &anc[0]
		root = NonCoding(root, substitutionMatrixFile, unitBranchLength)
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

	statePath := generateIlsStatePath(totalLength, numRoots, transMat, transMatConverted)
	statePathWig := stateToPostProbWig(statePath, numRoots, totalLength)
	topologyRecord := reconstruct.IlsWigToBed(statePathWig, chromName)
	topologyRecord, ilsEvolvedSeqs := reconstruct.SumIlsSeqsFasta(statePathWig, forwardEvolvedSeqs)

	return anc, forwardEvolvedSeqs, topologyRecord, ilsEvolvedSeqs
}

func pickState(prevState int, n int, transMatConverted mat.Matrix) int {
	prob := rand.Float64()
	for currentState := range n {
		if prob <= transMatConverted.At(prevState, currentState) {
			return currentState
		}
	}

	log.Fatal("Failed to pick next state.")
}

// pickInitialState samples an initial topology from the stationary distribution to avoid bias from a uniformly selected initial state.
func pickInitialState(transMat mat.Matrix) int {
	probs := matrix.StationaryDistribution(transMat)

	r := rand.Float64()
	cumulative := 0.0

	for i, prob := range probs {
		cumulative += prob
		if r < cumulative {
			return i
		}
	}

	return len(probs) - 1
}

// generateIlsStatePath simulates a sequence of topology states using a
// first-order Markov chain. The initial state is sampled from the stationary
// distribution of transMat, and subsequent states are sampled using
// transMatConverted, the cumulative transition matrix.
func generateIlsStatePath(totalLength int, n int, transMat mat.Matrix, transMatConverted mat.Matrix) []int {
	statePath := make([]int, totalLength)

	prevState := pickInitialState(transMat)
	statePath[0] = prevState

	for seqIdx := 1; seqIdx < totalLength; seqIdx++ {
		currentState := pickState(prevState, n, transMatConverted)
		statePath[seqIdx] = currentState
		prevState = currentState
	}

	return statePath
}

func stateToPostProbWig(states []int, numRoots int, seqLength int) map[string]wig.Wig {
	var topoName string
	sizes := make(map[string]chromInfo.ChromInfo, numRoots)
	for topo := range numRoots {
		topoName = fmt.Sprintf("V%d", topo)
		// add topo to map
		_, _ = sizes[topoName]
		sizes[topoName] = chromInfo.ChromInfo{Name: topoName, Size: seqLength, Order: topo}
	}

	stateWigs := wig.MakeSkeleton(sizes, 0) // return map[string]Wig
	for pos, state := range states {
		for topo := range numRoots {
			topoName = fmt.Sprintf("V%d", topo)
			if topo == state {
				stateWigs[topoName].Values[pos] = 1
			}
		}
	}

	return stateWigs
}
