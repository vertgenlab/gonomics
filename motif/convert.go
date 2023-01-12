package motif

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"math"
	"math/rand"
)

// PfmToPpm creates a position probability matrix from an input position frequency matrix.
// Pseudocounts may be applied for Laplace smoothing. The input float represents the value added
// to each cell of the Pfm. More info on pseudocounts; https://doi.org/10.1093/nar/gkn1019
func PfmToPpm(input PositionMatrix, pseudocount float64) PositionMatrix {
	if input.Type != Frequency {
		log.Fatalf("Input PositionMatrix must be of type 'Frequency' to be converted to a PPM.")
	}
	var row, column int
	var columnSum float64
	var answer = CopyPositionMatrix(input)
	answer.Type = Probability

	//for every column of the motif matrix
	for column = 0; column < len(input.Mat[0]); column++ {
		columnSum = input.Mat[0][column] + input.Mat[1][column] + input.Mat[2][column] + input.Mat[3][column] + (pseudocount * 4)
		for row = 0; row < 4; row++ { //for every row of that column of the motif matrix
			answer.Mat[row][column] = (input.Mat[row][column] + pseudocount) / columnSum
		}
	}
	return answer
}

// PfmSliceToPpmSlice creates a slice of position probabilty matrices
// from a slice of position frequency matrices.
func PfmSliceToPpmSlice(input []PositionMatrix, pseudocount float64) []PositionMatrix {
	var answer = make([]PositionMatrix, len(input))
	for i := range input {
		answer[i] = PfmToPpm(input[i], pseudocount)
	}
	return answer
}

// PpmToPwm creates a position weight matrix from an input position probability matrix.
func PpmToPwm(input PositionMatrix) PositionMatrix {
	if input.Type != Probability {
		log.Fatalf("Input PositionMatrix must be of type 'Probability' to be converted to a PWM.")
	}
	var row, column int
	var answer = CopyPositionMatrix(input)
	answer.Type = Weight
	for column = 0; column < len(input.Mat[0]); column++ {
		for row = 0; row < 4; row++ {
			answer.Mat[row][column] = math.Log2(input.Mat[row][column] * 4)
		}
	}
	return answer
}

// PpmSliceToPwmSlice creates a slice of position weight matrices from
// a slice of position probability matrices.
func PpmSliceToPwmSlice(input []PositionMatrix) []PositionMatrix {
	var answer = make([]PositionMatrix, len(input))
	for i := range input {
		answer[i] = PpmToPwm(input[i])
	}
	return answer
}

// ConsensusSequence takes in a PositionMatrix and returns the consensus motif sequence as a fasta.
// works for PFM, PPM, and PWM, as for all three, the consensus base is represented by the
// max column.
// if tieBreak is true, a tie in the PositionMatrix will produce a random output tied base.
func ConsensusSequence(input PositionMatrix, tieBreak bool) fasta.Fasta {
	var answer = make([]dna.Base, len(input.Mat[0]))
	var currMax int
	var currValue float64
	var row, column int
	for column = range input.Mat[0] {
		currMax = 0
		currValue = input.Mat[0][column]
		for row = 1; row < 4; row++ {
			if input.Mat[row][column] > currValue {
				currMax = row
				currValue = input.Mat[row][column]
			} else if tieBreak && input.Mat[row][column] == currValue && rand.Float64() > 0.5 { //in the case of a tie, we have a 50% chace of replacing the consensus output with the current base.
				currMax = row
			}
		}
		switch currMax {
		case 0:
			answer[column] = dna.A
		case 1:
			answer[column] = dna.C
		case 2:
			answer[column] = dna.G
		case 3:
			answer[column] = dna.T
		default:
			log.Fatalf("Error in ConsensusSequence. CurrMax value not recotnized.")
		}
	}
	return fasta.Fasta{Name: input.Name, Seq: answer}
}

// ConsensusSequences converts a slice of PositionMatrix structs into a slice of Fasta structs representing the consensus sequences for each matrix.
func ConsensusSequences(input []PositionMatrix, tieBreak bool) []fasta.Fasta {
	var answer = make([]fasta.Fasta, len(input))
	for i := range input {
		answer[i] = ConsensusSequence(input[i], tieBreak)
	}
	return answer
}

// ReverseComplement creates a new PositionMatrix representing the reverse complement
// position matrix of the input matrix.
func ReverseComplement(input PositionMatrix) PositionMatrix {
	//start by making a memory copy of the input array
	var answer = PositionMatrix{Id: input.Id, Name: input.Name, Mat: make([][]float64, 4)}
	for i := 0; i < 4; i++ {
		answer.Mat[i] = make([]float64, len(input.Mat[i]))
		copy(answer.Mat[i], input.Mat[i])
	}

	//now reverse the columns in every row
	var m, n int
	for row := 0; row < 4; row++ {
		for m, n = 0, len(answer.Mat[row])-1; m < n; m, n = m+1, n-1 {
			answer.Mat[row][m], answer.Mat[row][n] = answer.Mat[row][n], answer.Mat[row][m]
		}
	}

	//now complement the matrix
	answer.Mat[0], answer.Mat[3] = answer.Mat[3], answer.Mat[0]
	answer.Mat[1], answer.Mat[2] = answer.Mat[2], answer.Mat[1]

	return answer
}

// ReverseComplementAll creates a new slice of PositionMatrix structs, where each entry
// is the reverse complement position matrix of the corresponding index of the input slice
// of PositionMatrix structs.
func ReverseComplementAll(input []PositionMatrix) []PositionMatrix {
	var answer = make([]PositionMatrix, len(input))
	for i := range input {
		answer[i] = ReverseComplement(input[i])
	}
	return answer
}

// CopyPositionMatrix provides a memory copy of an input PositionMatrix struct.
func CopyPositionMatrix(input PositionMatrix) PositionMatrix {
	var row, column int
	var answer = PositionMatrix{Id: input.Id, Name: input.Name, Type: Probability, Mat: make([][]float64, 4)}
	for row = 0; row < 4; row++ {
		answer.Mat[row] = make([]float64, len(input.Mat[0]))
		for column = range answer.Mat[row] {
			answer.Mat[row][column] = input.Mat[row][column]
		}
	}

	return answer
}
