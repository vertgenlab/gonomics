package simulate

import (
	"fmt"
	"log"
	"math/rand"
	"strconv"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
)

// Main function to be called in simulateWrightFisher
func SimulateWrightFisher(set popgen.WrightFisherSettings) popgen.WrightFisherPopData {
	// Check if rFitness value is valid
	if set.RFitness < 0 {
		log.Fatalf("rFitness value must be greater or equal than zero. Found: %v.", set.RFitness)
	}

	if set.Verbose {
		fmt.Printf(`Population Size = %v
Genome Size = %v
Number of Generations = %v
Mutation Rate = %v
Relative Fitness = %v
GC Content = %v`, set.PopSize, set.GenomeSize, set.NumGen, set.MutRate, set.RFitness, set.GcContent)
	}
	// Randomly generate the initial sequence and copy to all individual in 2 fasta slices: curFasta and nextFasta
	curFasta, nextFasta := makeInitialPop(set)

	if set.Verbose {
		fmt.Printf("\nInitial sequence is\n%v\n", curFasta[0].Seq)
	}

	// Generate predetermined/stochastic relative fitness array based on fitness input value
	relFitArray := makeFitnessArray(curFasta[0].Seq, set)

	if set.Verbose {
		fmt.Printf("Relative fitness landscape is\n%v\n", relFitArray)
	}

	// An array keeping track of allele frequencies of all possible allelic states of all sites of all generation
	allFreq := makeAlleleFreqArray(curFasta, set)

	// Store the ancestral allele for each site
	ancestralAlleles := makeAncestralArray(curFasta[0].Seq, set)

	simulateAllGeneration(curFasta, nextFasta, relFitArray, allFreq, set)

	if set.Verbose {
		fmt.Printf("Ancestral alleles are\n%v\n", ancestralAlleles)
		fmt.Printf("Allele frequency landscape is\n%v\n", allFreq)
		fmt.Printf("Output Fasta is\n%v\n", curFasta)
	}

	wf := popgen.WrightFisherPopData{
		Fasta:     curFasta,
		Freq:      allFreq,
		Settings:  set,
		Ancestral: ancestralAlleles,
	}

	wf.Meta = makeMetadata(set)

	return wf
}

func makeInitialPop(set popgen.WrightFisherSettings) ([]fasta.Fasta, []fasta.Fasta) {
	curFasta := make([]fasta.Fasta, set.PopSize)
	nextFasta := make([]fasta.Fasta, set.PopSize)
	initialSeq := RandIntergenicSeq(set.GcContent, set.GenomeSize)

	for i := 0; i < set.PopSize; i++ {
		curFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
		curFasta[i].Seq = make([]dna.Base, set.GenomeSize)
		copy(curFasta[i].Seq, initialSeq)

		nextFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
		nextFasta[i].Seq = make([]dna.Base, set.GenomeSize)
		copy(nextFasta[i].Seq, initialSeq)
	}

	return curFasta, nextFasta
}

/*
 Helper function that constructs a slice that stores relative fitness of derived alleles on each position
 Think of it as a 2D array where row = position, column index = 0:A, 1:C, 2:G, 3:T
 However, I implemented only one allocation for the slice with size nrow * ncol
 And to access relative fitness of C allele at position 3, call slice[4*3 + 1]
*/
func makeFitnessArray(initSeq []dna.Base, set popgen.WrightFisherSettings) [][]float64 {
	answer := make([][]float64, set.GenomeSize)
	bases := [4]dna.Base{dna.A, dna.C, dna.G, dna.T}
	var i, j int
	//var r float64

	for i = 0; i < set.GenomeSize; i++ {
		answer[i] = make([]float64, 4)
		for j = 0; j < 4; j++ {
			if bases[j] == initSeq[i] { // check if this is ancestral allele
				answer[i][j] = float64(1) // relative fitness of ancestral allele is always 1
			} else {
				//r = rand.NormFloat64()*((rFitness-1)/3) + rFitness // this makes it that 3sd (99.7%) of generated fitness falls between the rFitness and 1
				//answer[i*4+j] = r
				answer[i][j] = set.RFitness
			}
		}
	}
	return answer
}

func makeAlleleFreqArray(curFasta []fasta.Fasta, set popgen.WrightFisherSettings) [][][]float64 {
	answer := make([][][]float64, set.NumGen+1)
	var t, s int
	for t = 0; t <= set.NumGen; t++ {
		answer[t] = make([][]float64, set.GenomeSize)
		for s = 0; s < set.GenomeSize; s++ {
			answer[t][s] = make([]float64, 4)
		}
		updateFreqArray(curFasta, t, answer)
	}
	return answer
}

func makeAncestralArray(initSeq []dna.Base, set popgen.WrightFisherSettings) []string {
	answer := make([]string, set.GenomeSize)

	for i := 0; i < set.GenomeSize; i++ {
		answer[i] = dna.BaseToString(initSeq[i])
	}
	return answer
}

func simulateAllGeneration(curFasta []fasta.Fasta, nextFasta []fasta.Fasta, relFitArray [][]float64, allFreq [][][]float64, set popgen.WrightFisherSettings) {
	var t, s, b, p int
	var r float64
	samplingPQRS := make([]float64, 4)
	normFactorArray := make([]float64, set.GenomeSize)
	updateNormFactorArray(0, relFitArray, allFreq, normFactorArray)

	// The below for loop could be written in another helper function to make this function more readable
	// 1st loop through every generation
	for t = 1; t <= set.NumGen; t++ {

		// 2nd loop through every base
		for s = 0; s < set.GenomeSize; s++ {

			// Calculate the allele frequency post-selection based on original frequency and relative fitness
			// sumFreqFitArray is the denominator that normalize the weighted frequency, rendering the sum of new frequencies = 1
			for b = 0; b < 4; b++ {
				samplingPQRS[b] = (allFreq[t-1][s][b] * relFitArray[s][b]) / (normFactorArray[s])
			}

			// fmt.Printf("Before Sampling = \n%v\nAfter Sampling = \n%v\n", allFreq[t-1][s], samplingPQRS)

			// 3rd loop through every individual
			for p = 0; p < set.PopSize; p++ {
				// Which chromosome from previous generation that pth individual inherits from?
				r = rand.Float64()
				if r < samplingPQRS[0] {
					nextFasta[p].Seq[s] = dna.A
				} else if r < sumSlice(samplingPQRS[0:2]) {
					nextFasta[p].Seq[s] = dna.C
				} else if r < sumSlice(samplingPQRS[0:3]) {
					nextFasta[p].Seq[s] = dna.G
				} else {
					nextFasta[p].Seq[s] = dna.T
				}

				// If random float generated is less than mutation rate (mutR), mutate the base
				if rand.Float64() < set.MutRate {
					nextFasta[p].Seq[s] = changeBase(nextFasta[p].Seq[s])
				}

			}
		}

		// Update the curFasta for the next generation to be nextFasta of this generation
		curFasta, nextFasta = nextFasta, curFasta
		updateFreqArray(curFasta, t, allFreq)
		updateNormFactorArray(t, relFitArray, allFreq, normFactorArray)
	}
}

func updateFreqArray(curFasta []fasta.Fasta, gen int, allFreq [][][]float64) {
	var a, c, g, t float64
	popSize := float64(len(curFasta))

	for s := 0; s < len(allFreq[gen]); s++ {
		a, c, g, t = 0, 0, 0, 0
		for i := 0; i < len(curFasta); i++ {
			switch curFasta[i].Seq[s] {
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
		allFreq[gen][s][0] = a / popSize
		allFreq[gen][s][1] = c / popSize
		allFreq[gen][s][2] = g / popSize
		allFreq[gen][s][3] = t / popSize
	}
}

/*
This helper function updates the sum of weighted relative fitnesset. This value is used for updating new allele frequencies after selection.
*/
func updateNormFactorArray(gen int, relFitArray [][]float64, allFreq [][][]float64, normFactorArray []float64) {
	for s := 0; s < len(normFactorArray); s++ {
		normFactorArray[s] = allFreq[gen][s][0]*relFitArray[s][0] +
			allFreq[gen][s][1]*relFitArray[s][1] +
			allFreq[gen][s][2]*relFitArray[s][2] +
			allFreq[gen][s][3]*relFitArray[s][3]
	}
}

func makeMetadata(set popgen.WrightFisherSettings) []string {
	meta := []string{
		"##PopulationSize=" + fmt.Sprint(set.PopSize),
		"NumGeneration=" + fmt.Sprint(set.NumGen),
		"Replicates=" + fmt.Sprint(set.GenomeSize),
		"MutationRate=" + strconv.FormatFloat(set.MutRate, 'g', 3, 64),
		"RelativeFitness=" + strconv.FormatFloat(set.RFitness, 'f', 5, 64),
	}

	return meta
}

// Helper function that sums the elements in slice
func sumSlice(slice []float64) float64 {
	var answer float64
	for _, v := range slice {
		answer += v
	}
	return answer
}
