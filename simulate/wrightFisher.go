package simulate

import (
	"fmt"
	"log"
	"math/rand"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
)

// Main function to be called in simulateWrightFisher.gp.
func SimulateWrightFisher(set popgen.WrightFisherSettings) popgen.WrightFisherPopData {
	checkValidInput(set)                          // Check various inputs
	set.AncestralAllele = setAncestralAllele(set) // Set ancestral allele if given by input

	if set.Verbose {
		fmt.Printf(`Population Size = %v
Genome Size = %v
Number of Generations = %v
Mutation Rate = %v
Relative Fitness = %v
GC Content = %v`, set.PopSize, set.GenomeSize, set.NumGen, set.MutRate, set.RFitness, set.GcContent)
	}

	// Make two slices of fasta containing initial sequences. Default: random generated sequences
	// If initFreq or fitnessString is specified, all sequences start with the specified ancestrall allele
	curFasta, nextFasta := makeInitialPop(set)

	if set.Verbose {
		fmt.Printf("\nInitial sequence is\n%v\n", curFasta[0].Seq)
	}

	// An array keeping track of allele frequencies of all possible allelic states of all sites of all generation
	allFreq := makeAlleleFreqArray(curFasta, set)

	// Store the ancestral allele for each site in the sequence
	ancestralAlleles := makeAncestralArray(curFasta[0].Seq, set)

	// Generate the relative fitness array based on fitness input value
	relFitArray := makeFitnessArray(curFasta[0].Seq, set)

	if set.Verbose {
		fmt.Printf("Relative fitness landscape is\n%v\n", relFitArray)
	}

	// Main function that loops through all generation, individual, and site
	simulateAllGeneration(curFasta, nextFasta, relFitArray, allFreq, set)

	if set.Verbose {
		fmt.Printf("Ancestral alleles are\n%v\n", ancestralAlleles)
		fmt.Printf("Allele frequency landscape is\n%v\n", allFreq)
		fmt.Printf("Output Fasta is\n%v\n", curFasta)
	}

	// The output is a WrightFisherPopData struct (see popgen)
	wf := popgen.WrightFisherPopData{
		Fasta:     curFasta,
		Freq:      allFreq,
		Settings:  set,
		Ancestral: ancestralAlleles,
	}

	// Generate metadata for the output file
	wf.Meta = makeMetadata(set)

	return wf
}

/*
checkValidInput() checks if there are inputs from flags that are invalid (value or format), if yes, stop the program.
*/
func checkValidInput(set popgen.WrightFisherSettings) {
	if set.InitFreq != "" && set.FitnessString != "" {
		fAncestral := strings.ToUpper(strings.Split(set.FitnessString, ",")[4])
		iAncestral := strings.ToUpper(strings.Split(set.InitFreq, ",")[4])
		if fAncestral != iAncestral {
			log.Fatalf("Ancestral alleles in -i and -W must be the same: %v != %v.", iAncestral, fAncestral)
		}
	}
	if set.RFitness < 0 {
		log.Fatalf("rFitness value must be greater or equal than zero. Found: %v.", set.RFitness)
	}
}

/*
setAncestralAllele() returns what the one ancestral allele is if given by input.
*/
func setAncestralAllele(set popgen.WrightFisherSettings) string {
	var answer string
	if set.InitFreq == "" && set.FitnessString == "" {
		answer = ""
	} else if set.InitFreq != "" {
		answer = strings.ToUpper(strings.Split(set.InitFreq, ",")[4])
	} else {
		answer = strings.ToUpper(strings.Split(set.FitnessString, ",")[4])
	}
	return answer
}

/*
makeInitialPop() returns identical slices of fasta containing initial sequences
Case I-no ancestral allele is specified: randomly generate sequences
Case II-ancestral allele is specified: generate repeats of that allele.
*/
func makeInitialPop(set popgen.WrightFisherSettings) ([]fasta.Fasta, []fasta.Fasta) {
	curFasta := make([]fasta.Fasta, set.PopSize)
	nextFasta := make([]fasta.Fasta, set.PopSize)
	var i int

	// Case I: No ancestral allele is specified. Random generate sequence.
	if set.AncestralAllele == "" {
		initialSeq := RandIntergenicSeq(set.GcContent, set.GenomeSize)
		for i = 0; i < set.PopSize; i++ {
			curFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
			curFasta[i].Seq = make([]dna.Base, set.GenomeSize)
			copy(curFasta[i].Seq, initialSeq)

			nextFasta[i].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(i))
			nextFasta[i].Seq = make([]dna.Base, set.GenomeSize)
			copy(nextFasta[i].Seq, initialSeq)
		}
	} else { // Case II: Ancestral allele is specified:
		a := strings.Split(set.InitFreq, ",")
		freq := make([]float64, 4)
		var e error

		// Parse all string input of frequency into float64
		for i = 0; i < 4; i++ {
			freq[i], e = strconv.ParseFloat(a[i], 64)
			if e != nil {
				fmt.Println("Invalid initial frequencies")
				exception.PanicOnErr(e)
			}
		}

		// Fatal if the sum of frequencies is not equal to 1
		if sumSlice(freq) != 1.0 {
			log.Fatalf("The sum of initial frequencies must be 1")
		}

		// Make fasta based on the given allele frequencies
		makeFastaByRatio(curFasta, nextFasta, freq, set)

	}
	return curFasta, nextFasta
}

/*
makeFastaByRatio() makes cur- and nextFasta from the given input frequencies.
*/
func makeFastaByRatio(curFasta []fasta.Fasta, nextFasta []fasta.Fasta, freq []float64, set popgen.WrightFisherSettings) {
	var ratio float64 // pointer keeps track of current ratio of population that is filled by an allele
	for j := 0; j < set.PopSize; j++ {
		curFasta[j].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(j))
		curFasta[j].Seq = make([]dna.Base, set.GenomeSize)

		nextFasta[j].Name = fmt.Sprintf("Seq_%v", strconv.Itoa(j))
		nextFasta[j].Seq = make([]dna.Base, set.GenomeSize)

		ratio = float64(j+1) / float64(set.PopSize)

		if ratio <= freq[0] {
			makeUniformSeq(curFasta, j, ratio, dna.A, set)
		} else if ratio <= sumSlice(freq[0:2]) {
			makeUniformSeq(curFasta, j, ratio, dna.C, set)
		} else if ratio <= sumSlice(freq[0:3]) {
			makeUniformSeq(curFasta, j, ratio, dna.G, set)
		} else {
			makeUniformSeq(curFasta, j, ratio, dna.T, set)
		}
		copy(nextFasta[j].Seq, curFasta[j].Seq)
	}
}

/*
makeUniformSeq() makes a sequence of uniform repeats of one allele.
*/
func makeUniformSeq(curFasta []fasta.Fasta, j int, ratio float64, base dna.Base, set popgen.WrightFisherSettings) {
	for k := 0; k < set.GenomeSize; k++ {
		curFasta[j].Seq[k] = base
	}
}

/*
makeAlleleFreqArray() makes an allele frequency array
3D array, zero-based, [generation][site][base].
*/
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

/*
makeAncestralArray() makes a slice that stores the ancestral alleles of every site.
*/
func makeAncestralArray(initSeq []dna.Base, set popgen.WrightFisherSettings) []string {
	answer := make([]string, set.GenomeSize)

	// Case I: Ancestral allele is not specified
	if set.InitFreq == "" {
		for i := 0; i < set.GenomeSize; i++ {
			answer[i] = dna.BaseToString(initSeq[i])
		}
	} else { // Case II: Ancecstral allele is specifie
		for i := 0; i < set.GenomeSize; i++ {
			answer[i] = set.AncestralAllele
		}
	}

	return answer
}

/*
makeFitnessArray() constructs a slice that stores relative fitness of derived alleles on each position
2D array, zero-based, row = site, column = 0:A, 1:C, 2:G, 3:T.
*/
func makeFitnessArray(initSeq []dna.Base, set popgen.WrightFisherSettings) [][]float64 {
	answer := make([][]float64, set.GenomeSize)
	bases := [4]dna.Base{dna.A, dna.C, dna.G, dna.T}
	var i, j int

	if set.FitnessString == "" {
		for i = 0; i < set.GenomeSize; i++ {
			answer[i] = make([]float64, 4)
			for j = 0; j < 4; j++ {
				if bases[j] == initSeq[i] { // check if this is ancestral allele
					answer[i][j] = float64(1) // relative fitness of ancestral allele is always 1
				} else {
					answer[i][j] = set.RFitness
				}
			}
		}
	} else {
		a := strings.Split(set.FitnessString, ",")
		fitness := make([]float64, 4)
		var e error

		// Parse all string input of fitness value into float64
		for i = 0; i < 4; i++ {
			fitness[i], e = strconv.ParseFloat(a[i], 64)
			if fitness[i] < 0 {
				log.Fatalf("Relative fitness values must be greater or equal than zero. Found: %v.", set.FitnessString)
			}
			if e != nil {
				fmt.Println("Invalid relative fitness string")
				exception.PanicOnErr(e)
			}

		}

		// Make the fitness array
		for i = 0; i < set.GenomeSize; i++ {
			answer[i] = make([]float64, 4)
			for j = 0; j < 4; j++ {
				answer[i][j] = fitness[j]
			}
		}
	}

	return answer
}

/*
simulateAllGeneration() simulates the changes in allele frequencies through all generation, all individual, and all site.
*/
func simulateAllGeneration(curFasta []fasta.Fasta, nextFasta []fasta.Fasta, relFitArray [][]float64, allFreq [][][]float64, set popgen.WrightFisherSettings) {
	var t, s, b, p int
	var r float64
	// This slice contains new frequencies of each allele after they are resampled based on relative fitness
	samplingPQRS := make([]float64, 4)
	normFactorArray := make([]float64, set.GenomeSize) // Normalizing factor for the new frequencies (sum = 1)
	// Set the initial normalizing factor (generation 0)
	updateNormFactorArray(0, relFitArray, allFreq, normFactorArray)

	// 1st loop through every generation
	for t = 1; t <= set.NumGen; t++ {

		// 2nd loop through every base
		for s = 0; s < set.GenomeSize; s++ {

			// Calculate the allele frequency post-selection based on original frequency and relative fitness
			// normFactorArray stores the denominator that normalize the weighted frequency, rendering the sum of new frequencies = 1
			for b = 0; b < 4; b++ {
				samplingPQRS[b] = (allFreq[t-1][s][b] * relFitArray[s][b]) / (normFactorArray[s])
			}

			// 3rd loop through every individual
			for p = 0; p < set.PopSize; p++ {
				// Determines which chromosome from previous generation that pth individual inherits from?
				// The probability is based on the frequency of the corresponding allele.
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

				// If random float generated is less than the set mutation rate, mutate the base (random)
				if rand.Float64() < set.MutRate {
					nextFasta[p].Seq[s] = mutate(nextFasta[p].Seq[s], set.GcContent)
				}
			}
		}

		// Update the curFasta for the next generation to be nextFasta of this generation
		curFasta, nextFasta = nextFasta, curFasta
		updateFreqArray(curFasta, t, allFreq)
		updateNormFactorArray(t, relFitArray, allFreq, normFactorArray)
	}
}

/*
updateFreqArray() adds the new set of frequencies into the allFreq based on the current slice of fasta for a specified generation
This is needed because we only have two slices of fasta to keep track of the sequences at the time (fasta has no memory of sequence >2 generations prior)
This update is called once every generation (because of independent sites, and we keep separate curFasta and nextFasta).
*/
func updateFreqArray(curFasta []fasta.Fasta, gen int, allFreq [][][]float64) {
	var a, c, g, t float64 // Counters for each base (float64 because of later division)
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

		// Calculate the new frequencies for next generation
		allFreq[gen][s][0] = a / popSize
		allFreq[gen][s][1] = c / popSize
		allFreq[gen][s][2] = g / popSize
		allFreq[gen][s][3] = t / popSize
	}
}

/*
updateNormFactorArray() calculates and updates the normalizing factor
based on the original frequencies before sampling and relative fitness.
*/
func updateNormFactorArray(gen int, relFitArray [][]float64, allFreq [][][]float64, normFactorArray []float64) {
	for s := 0; s < len(normFactorArray); s++ {
		normFactorArray[s] = allFreq[gen][s][0]*relFitArray[s][0] +
			allFreq[gen][s][1]*relFitArray[s][1] +
			allFreq[gen][s][2]*relFitArray[s][2] +
			allFreq[gen][s][3]*relFitArray[s][3]
	}
}

/*
makeMetadata() takes values from WrightFisherSettings, and turns into a slice of string with important parameters.
This is to be printed as the comments on top of output tsv.
*/
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

/*
sumSlice() sums the values inside a []float64.
*/
func sumSlice(slice []float64) float64 {
	var answer float64
	for _, v := range slice {
		answer += v
	}
	return answer
}

/*
There exists the similar function in simulate.go, but that function doesn't allow
me to choose GC content (it was set to 0.41 by default and not mutable).
*/
func mutate(originalBase dna.Base, GC float64) dna.Base {
	newBase := chooseRandomBase(GC)

	for newBase == originalBase {
		newBase = chooseRandomBase(GC)
	}
	return newBase
}
