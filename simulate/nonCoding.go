package simulate

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math/rand"
	"strings"
)

// Jukes-Cantor matrix, medium distance.
var defaultSubstitutionMatrix [][]float64 = [][]float64{
	{0.91, 0.03, 0.03, 0.03},
	{0.03, 0.91, 0.03, 0.03},
	{0.03, 0.03, 0.91, 0.03},
	{0.03, 0.03, 0.03, 0.91},
}

func NonCoding(root *expandedTree.ETree, substitutionMatrixFile string, unitBranchLength float64) *expandedTree.ETree {
	unitMatrix := ParseSubstitutionMatrix(substitutionMatrixFile)
	expandedTree.PopulateSubstitutionMatrices(root, unitMatrix, unitBranchLength)
	recursiveEvolveSequence(root)
	return root
}

// recursiveEvolveSequence is a helper function of NonCoding which performs a molecular evolution simulation on
// the two children of each node in an input tree, recursively.
func recursiveEvolveSequence(node *expandedTree.ETree) {
	if node.Left != nil {
		node.Left.Fasta = substituteSequenceWithMatrix(node.Fasta, node.Left.SubstitutionMatrix, node.Left.Name)
		recursiveEvolveSequence(node.Left)
	}
	if node.Right != nil {
		node.Right.Fasta = substituteSequenceWithMatrix(node.Fasta, node.Right.SubstitutionMatrix, node.Right.Name)
		recursiveEvolveSequence(node.Right)
	}
}

// substituteSequenceWithMatrix runs a forward substitution simulation on an input fasta sequence, based on
// a defined 4x4 substitution matrix.
func substituteSequenceWithMatrix(inSeq *fasta.Fasta, substitutionMatrix [][]float64, outSeqName string) *fasta.Fasta {
	var answer fasta.Fasta = fasta.Fasta{Name: outSeqName, Seq: make([]dna.Base, len(inSeq.Seq))}
	for currPos := range inSeq.Seq {
		substitutedBase := substituteWithMatrix(inSeq.Seq[currPos], substitutionMatrix)
		answer.Seq[currPos] = substitutedBase
	}
	return &answer
}

// substitutionWithMatrix passes an input dna.Base through a substitutionMatrix. This simulates the mutation process
// for a given substitutionMatrix for a single input base.
func substituteWithMatrix(inBase dna.Base, substitutionMatrix [][]float64) dna.Base {
	var currRand float64 = rand.Float64()
	if inBase > 3 {
		return inBase
	}
	if currRand < substitutionMatrix[inBase][dna.A] {
		return dna.A
	} else if currRand < substitutionMatrix[inBase][dna.A]+substitutionMatrix[inBase][dna.C] {
		return dna.C
	} else if currRand < substitutionMatrix[inBase][dna.A]+substitutionMatrix[inBase][dna.C]+substitutionMatrix[inBase][dna.G] {
		return dna.G
	}
	return dna.T
}

// ParseSubstitutionMatrix reads a substitution matrix from an input file and returns it as a [][]float64
func ParseSubstitutionMatrix(filename string) [][]float64 {
	if filename == "" {
		return defaultSubstitutionMatrix
	}
	lines := fileio.Read(filename)
	if len(lines) != 4 {
		log.Fatalf("Error: expected 4 lines in input substitution matrix. Found: %v.\n", len(lines))
	}
	var answer = make([][]float64, 4)
	var words []string
	var currColumn int
	for currRow := range lines {
		answer[currRow] = make([]float64, 4)
		words = strings.Split(lines[currRow], "\t")
		if len(words) != 4 {
			log.Fatalf("Error: unable to parse substitution matrix file. Expected 4 fields per line:\n%v\n", lines[currRow])
		}
		for currColumn = range words {
			answer[currRow][currColumn] = parse.StringToFloat64(words[currColumn])
		}
	}
	return answer
}
