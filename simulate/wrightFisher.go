package simulate

import (
	"fmt"
	"math/rand"
	"strconv"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

// Main function to be called in simulateWrightFisher
func SimulateWrightFisher(popSize int, mutR float64, numGen int, genomeSize int, rFitness float64, gcContent float64, verbose int) []fasta.Fasta {
	curFasta := make([]fasta.Fasta, popSize)
	nextFasta := make([]fasta.Fasta, popSize)

	ancestralAlleles := make([]dna.Base, genomeSize)
	aFreqArray := make([]float64, genomeSize*4)

	sumFreqFitArray := make([]float64, genomeSize)

	samplingPQRS := make([]float64, 4) // Array containing allele frequencies for each allele at sampling (after selection)

	//fmt.Println(byte(65) == 'A')
	//fmt.Println(dna.ByteToBase(uint8(0)))
	var i, t, b, p int
	var r float64

	// Randomly generate the initial sequence and copy to all individual in both curFasta and nextFasta
	makeInitialPop(curFasta, nextFasta, popSize, genomeSize, gcContent)

	relFitArray := makeFitnessArray(genomeSize, rFitness, curFasta[0].Seq)
	copy(ancestralAlleles, curFasta[0].Seq)

	for b = 0; b < genomeSize; b++ {
		updateAlleleFreqArray(curFasta, b, aFreqArray)
		updateSumFreqFitArray(b, relFitArray, aFreqArray, sumFreqFitArray)
	}

	// 1st loop through every generation
	for t = 1; t < numGen; t++ {

		// 2nd loop through every base
		for b = 0; b < genomeSize; b++ {

			// Calculate the allele frequency post-selection based on original frequency and relative fitness
			// sumFreqFitArray is the denominator that normalize the weighted frequency, rendering the sum of new frequencies = 1
			for i = 0; i < 4; i++ {
				samplingPQRS[i] = (aFreqArray[b*4+i] * relFitArray[b*4+1]) / (sumFreqFitArray[b])
			}

			// 3rd loop through every individual
			for p = 0; p < popSize; p++ {
				// Which chromosome from previous generation that pth individual inherits from?

				r = rand.Float64()
				if r < samplingPQRS[0] {
					nextFasta[p].Seq[b] = dna.A
				} else if r < sumSlice(samplingPQRS[0:2]) {
					nextFasta[p].Seq[b] = dna.C
				} else if r < sumSlice(samplingPQRS[0:3]) {
					nextFasta[p].Seq[b] = dna.G
				} else {
					nextFasta[p].Seq[b] = dna.T
				}

				// If random float generated is less than mutation rate (mutR), mutate the base
				if rand.Float64() < mutR {
					nextFasta[p].Seq[b] = changeBase(nextFasta[p].Seq[b])
				}

			}
			// This helper function only updates one position at a time
			updateAlleleFreqArray(nextFasta, b, aFreqArray)
			updateSumFreqFitArray(b, relFitArray, aFreqArray, sumFreqFitArray)
		}

		// Update the curFasta for the next generation to be nextFasta of this generation
		curFasta, nextFasta = nextFasta, curFasta
	}

	if verbose == 1 {
		fmt.Println("Allele Frequency")
		fmt.Println(aFreqArray)
		fmt.Println("Relative Fitness Array")
		fmt.Println(relFitArray)
		fmt.Println("Sum Frequency * Fitness")
		fmt.Println(sumFreqFitArray)
		fmt.Println("Output Fasta")
		fmt.Println(curFasta)
	}

	return curFasta
}

func makeInitialPop(curFasta []fasta.Fasta, nextFasta []fasta.Fasta, popSize int, genomeSize int, gcContent float64) {
	var i int

	initialSeq := RandIntergenicSeq(gcContent, genomeSize)

	for i = 0; i < popSize; i++ {
		curFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
		curFasta[i].Seq = make([]dna.Base, genomeSize)
		copy(curFasta[i].Seq, initialSeq)

		nextFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
		nextFasta[i].Seq = make([]dna.Base, genomeSize)
		copy(nextFasta[i].Seq, initialSeq)
	}
}

/*
 Helper function that constructs a slice that stores relative fitness of derived alleles on each position
 Think of it as a 2D array where row = position, column index = 0:A, 1:C, 2:G, 3:T
 However, I implemented only one allocation for the slice with size nrow * ncol
 And to access relative fitness of C allele at position 3, call slice[4*3 + 1]
*/
func makeFitnessArray(genomeSize int, rFitness float64, initSeq []dna.Base) []float64 {
	answer := make([]float64, genomeSize*4)
	var i, j int
	//var r float64

	bases := [4]dna.Base{dna.A, dna.C, dna.G, dna.T}

	for i = 0; i < genomeSize; i++ {
		for j = 0; j < 4; j++ {
			if bases[j] == initSeq[i] { // check if this is ancestral allele
				answer[i*4+j] = float64(1) // relative fitness of ancestral allele is always 1
			} else {
				//r = rand.NormFloat64()*((rFitness-1)/3) + rFitness // this makes it that 3sd (99.7%) of generated fitness falls between the rFitness and 1
				//answer[i*4+j] = r
				answer[i*4+j] = rFitness
			}
		}
	}
	return answer
}

/*
This helper function updates the allele frequency array AFTER all new generation is sampled and has gone through meiosis (whether to have mutation or not)
*/
func updateAlleleFreqArray(curFasta []fasta.Fasta, pos int, aFreqArray []float64) {
	var i int
	var a, c, g, t float64
	for i = 0; i < len(curFasta); i++ {
		switch curFasta[i].Seq[pos] {
		case dna.A:
			a++
		case dna.C:
			c++
		case dna.G:
			g++
		case dna.T:
			t++
		}
	}

	aFreqArray[pos*4+0] = a / float64(len(curFasta))
	aFreqArray[pos*4+1] = c / float64(len(curFasta))
	aFreqArray[pos*4+2] = g / float64(len(curFasta))
	aFreqArray[pos*4+3] = t / float64(len(curFasta))
}

func updateSumFreqFitArray(pos int, relFitArray []float64, aFreqArray []float64, sumFreqFitArray []float64) {
	sumFreqFitArray[pos] = aFreqArray[pos*4+0]*relFitArray[pos*4+0] +
		aFreqArray[pos*4+1]*relFitArray[pos*4+1] +
		aFreqArray[pos*4+2]*relFitArray[pos*4+2] +
		aFreqArray[pos*4+3]*relFitArray[pos*4+3]
}

// Helper function that sums the elements in slice
func sumSlice(slice []float64) float64 {
	var answer float64
	for _, v := range slice {
		answer += v
	}
	return answer
}
