package jaspar

// PfmToPpm converts an input position frequency matrix to a position probability matrix.
// Pseudocounts may be applied for laplace smoothing. The input float represents the value added
// to each cell of the Pfm. More info on pseudocounts; https://doi.org/10.1093/nar/gkn1019
func PfmToPpm(input Pfm, pseudocount float64) Ppm {
	var answer = Ppm{Id: input.Id, Name: input.Name, Mat: make([][]float64, 4)}
	var columnSum float64
	for i := 0; i < 4; i ++ {
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
func PfmSliceToPpmSlice(input []Pfm, pseudocount float64) []Ppm {
	var answer = make([]Ppm, 0)
	for i := range input {
		answer = append(answer, PfmToPpm(input[i], pseudocount))
	}
	return answer
}
