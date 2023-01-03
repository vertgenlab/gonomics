package motif

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"math"
	"math/rand"
)

// PfmToPpm converts an input position frequency matrix to a position probability matrix.
// Pseudocounts may be applied for laplace smoothing. The input float represents the value added
// to each cell of the Pfm. More info on pseudocounts; https://doi.org/10.1093/nar/gkn1019
func PfmToPpm(input PositionMatrix, pseudocount float64) PositionMatrix {
	if input.Type != Frequency {
		log.Fatalf("Input PositionMatrix must be of type 'Frequency' to be converted to a PPM.")
	}
	var answer = PositionMatrix{Id: input.Id, Name: input.Name, Type: Probability, Mat: make([][]float64, 4)}
	var columnSum float64
	for i := 0; i < 4; i++ {
		answer.Mat[i] = make([]float64, len(input.Mat[0]))
	}

	//for every column of the motif matrix
	for i := 0; i < len(input.Mat[0]); i++ {
		columnSum = input.Mat[0][i] + input.Mat[1][i] + input.Mat[2][i] + input.Mat[3][i] + (pseudocount * 4)
		for j := 0; j < 4; j++ { //for every row of that column of the motif matrix
			answer.Mat[j][i] = (input.Mat[j][i] + pseudocount) / columnSum
		}
	}
	return answer
}

// PfmSliceToPpmSlice converts a slice of position frequency matrices to
// a slice of position probability matrices.
func PfmSliceToPpmSlice(input []PositionMatrix, pseudocount float64) []PositionMatrix {
	var answer = make([]PositionMatrix, 0)
	for i := range input {
		answer = append(answer, PfmToPpm(input[i], pseudocount))
	}
	return answer
}

// PpmToPwm converts an input position probability matrix to a position weight matrix.
func PpmToPwm(input PositionMatrix) PositionMatrix {
	if input.Type != Probability {
		log.Fatalf("Input PositionMatrix must be of type 'Probability' to be converted to a PWM.")
	}
	var answer = PositionMatrix{Id: input.Id, Name: input.Name, Type: Weight, Mat: make([][]float64, 4)}
	for i := 0; i < 4; i++ {
		answer.Mat[i] = make([]float64, len(input.Mat[0]))
	}
	for i := 0; i < len(input.Mat[0]); i++ {
		for j := 0; j < 4; j++ {
			answer.Mat[j][i] = math.Log2(input.Mat[j][i] * 4)
		}
	}
	return answer
}

// PpmSliceToPwmSlice converts a slice of position probability matrices to
// a slice of position weight matrices.
func PpmSliceToPwmSlice(input []PositionMatrix) []PositionMatrix {
	var answer = make([]PositionMatrix, 0)
	for i := range input {
		answer = append(answer, PpmToPwm(input[i]))
	}
	return answer
}

// ConsensusSequence takes in a PositionMatrix and returns the consensus motif sequence as a fasta.
// works for PFM, PPM, and PWM, as for all three, the consensus base is represented by the
// max column.
// if tieBreak is true, a tie in the PositionMatrix will produce a random output tied base.
func ConsensusSequence(input PositionMatrix, tieBreak bool) fasta.Fasta {
	var answer = make([]dna.Base, 0)
	var currMax int
	var currValue float64
	for i := range input.Mat[0] {
		currMax = 0
		currValue = input.Mat[0][i]
		for j := 1; j < 4; j++ {
			if input.Mat[j][i] > currValue {
				currMax = j
				currValue = input.Mat[j][i]
			} else if tieBreak && input.Mat[j][i] == currValue && rand.Float64() > 0.5 { //in the case of a tie, we have a 50% chace of replacing the consensus output with the current base.
				currMax = j
			}
		}
		switch currMax {
		case 0:
			answer = append(answer, dna.A)
		case 1:
			answer = append(answer, dna.C)
		case 2:
			answer = append(answer, dna.G)
		case 3:
			answer = append(answer, dna.T)
		default:
			log.Fatalf("Error in ConsensusSequence. CurrMax value not recotnized.")
		}
	}
	return fasta.Fasta{Name: input.Name, Seq: answer}
}

// ConsensusSequences converts a slice of PositionMatrix structs into a slice of Fasta structs representing the consensus sequences for each matrix.
func ConsensusSequences(input []PositionMatrix, tieBreak bool) []fasta.Fasta {
	var answer = make([]fasta.Fasta, 0)
	for i := range input {
		answer = append(answer, ConsensusSequence(input[i], tieBreak))
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
	var answer = make([]PositionMatrix, 0)
	for i := range input {
		answer = append(answer, ReverseComplement(input[i]))
	}
	return answer
}
