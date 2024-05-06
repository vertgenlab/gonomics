package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

// PriorSettings defines the set of options and arguments for the samAssembler prior subcommand.
type PriorSettings struct {
	SamFileName         string
	ReferenceFile       string
	OutFile             string
	Epsilon             float64
	LikelihoodCacheSize int
	PseudoCount         float64
	MinCoverage         int
	AsCounts            bool
}

// EmpiricalErrorEstimator keeps track of several variables used to estimate sequencing error rates (epsilon)
// and cytosine deamination error rates (lambda) from pileup data.
type EmpiricalErrorEstimator struct {
	NumEpsilon             int
	NumLambdaPlusEpsilon   int
	TotalEpsilon           int
	TotalLambdaPlusEpsilon int
}

// priorUsage defines the help message for the samAssembler prior subcommand and prints default options.
func priorUsage(priorFlags *flag.FlagSet) {
	fmt.Print("samAssembler prior - Construct an empirical conditional Dirichlet prior for output diploid" +
		"genotypes based on maximum likelihood estimation of genotypes from aligned short reads.\n" +
		"Usage: \n" +
		"samAssembler [options] prior reads.sam/bam ref.fa output.txt\n" +
		"Options:\n")
	priorFlags.PrintDefaults()
}

// parsePriorArgs is the main function of the samAssembler prior subcommand. Defines and parses arguments before running the SamAssemblerPrior function.
func parsePriorArgs() {
	var err error
	var expectedNumArgs int = 3
	priorFlags := flag.NewFlagSet("prior", flag.ExitOnError)

	var epsilon *float64 = priorFlags.Float64("epsilon", 0.01, "Set the expected misclassification error rate.")
	var likelihoodCacheSize *int = priorFlags.Int("likelihoodCacheSize", 100, "Set the maximum dimension of the likelihood caches. Should be slightly larger than highest expected pile depth.")
	var pseudoCount *float64 = priorFlags.Float64("pseudoCount", 0.01, "Prime the empirical prior with pseudoCounts to avoid prior probabilities of 0.\n"+
		"Increasing this value regularizes the prior against overfitting.\n"+
		"The model will underfit, biased towards Jukes-Cantor, if this value is set too high.")
	var asCounts *bool = priorFlags.Bool("asCounts", false, "Return output as counts, instead of probabilities. Useful for debugging, but count matrices cannot be used directly in 'build'.")
	var minCoverage *int = priorFlags.Int("minCoverage", 0, "Specifies a threshold where only piles with coverage above this value contribute to the substitution matrix and error  rate estimations.")

	err = priorFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	priorFlags.Usage = func() { priorUsage(priorFlags) }

	if len(priorFlags.Args()) != expectedNumArgs {
		priorFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(priorFlags.Args()))
	}

	samFileName := priorFlags.Arg(0)
	referenceFileName := priorFlags.Arg(1)
	outFileName := priorFlags.Arg(2)

	s := PriorSettings{
		SamFileName:         samFileName,
		ReferenceFile:       referenceFileName,
		OutFile:             outFileName,
		Epsilon:             *epsilon,
		LikelihoodCacheSize: *likelihoodCacheSize,
		PseudoCount:         *pseudoCount,
		AsCounts:            *asCounts,
		MinCoverage:         *minCoverage,
	}
	SamAssemblerPrior(s)
}

// SamAssemblerPrior creates an empirical prior from SAM/BAM pileup data for the samAssembler build command.
func SamAssemblerPrior(s PriorSettings) {
	// substitution matrix for priors
	var answer = make([][]float64, 4)
	// set pseudocount in matrix
	for i := range answer {
		answer[i] = make([]float64, 10)
		for j := range answer[i] {
			answer[i][j] = s.PseudoCount
		}
	}

	var errorEstimator = EmpiricalErrorEstimator{NumEpsilon: 0, NumLambdaPlusEpsilon: 0, TotalEpsilon: 0, TotalLambdaPlusEpsilon: 0}

	var refBase dna.Base
	var currChrom string
	var baseCall sam.DiploidBase
	var err error

	// read pileups from sam/bam
	reads, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)
	//read reference from fasta
	ref := fasta.Read(s.ReferenceFile)
	for i := range ref {
		dna.AllToUpper(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)
	// cache of substitution matrix if assembling genome is homozygous
	homozygousCache := make([][]float64, s.LikelihoodCacheSize)
	for i := range homozygousCache {
		homozygousCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	// cache of substitution matrix if assembling genome is heterozygous
	heterozygousCache := make([][]float64, s.LikelihoodCacheSize)
	for i := range heterozygousCache {
		heterozygousCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	// non-empirical cache, just flat, equally probable substitutions between all bases
	diploidBasePriorCache := sam.MakeDiploidBaseFlatPriorCache()
	// for all pileups
	for p := range piles {
		currChrom = header.Chroms[p.RefIdx].Name
		refBase = refMap[currChrom][p.Pos-1]
		if getPileCoverage(p) > s.MinCoverage {
			// if the baseCall doesn't return a genotype with an unknown base, then add the base to the transition matrix
			if refBase < 4 {
				baseCall = sam.DiploidBaseCallFromPile(p, refBase, diploidBasePriorCache, homozygousCache, heterozygousCache, sam.AncientLikelihoodCache{}, s.Epsilon, 0)
				if baseCall < 10 { //checks that we have a valid genotype (not NN)
					answer[refBase][baseCall]++
					updateErrorEstimate(&errorEstimator, baseCall, p)
				}
			}
		}
	}

	// epsilon MLE is just the empirical proportion of reads in homozygous sites that have incorrect bases.
	EpsilonEstimate := float64(errorEstimator.NumEpsilon) / float64(errorEstimator.TotalEpsilon)
	LambdaEstimate := numbers.Max((float64(errorEstimator.NumLambdaPlusEpsilon)/float64(errorEstimator.TotalLambdaPlusEpsilon))-EpsilonEstimate, 0)

	// convert matrix to probability
	if !s.AsCounts {
		answer = convertToProb(answer)
	}
	// format output file using aliases
	out := fileio.EasyCreate(s.OutFile)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Epsilon\t%v\n", EpsilonEstimate)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Lambda\t%v\n", LambdaEstimate)
	_, err = fmt.Fprintf(out, ".\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n")
	// 10 possible genotypes per reference genotype
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefA\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.A][0],
		answer[dna.A][1],
		answer[dna.A][2],
		answer[dna.A][3],
		answer[dna.A][4],
		answer[dna.A][5],
		answer[dna.A][6],
		answer[dna.A][7],
		answer[dna.A][8],
		answer[dna.A][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefC\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.C][0],
		answer[dna.C][1],
		answer[dna.C][2],
		answer[dna.C][3],
		answer[dna.C][4],
		answer[dna.C][5],
		answer[dna.C][6],
		answer[dna.C][7],
		answer[dna.C][8],
		answer[dna.C][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefG\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.G][0],
		answer[dna.G][1],
		answer[dna.G][2],
		answer[dna.G][3],
		answer[dna.G][4],
		answer[dna.G][5],
		answer[dna.G][6],
		answer[dna.G][7],
		answer[dna.G][8],
		answer[dna.G][9],
	)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "RefT\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n",
		answer[dna.T][0],
		answer[dna.T][1],
		answer[dna.T][2],
		answer[dna.T][3],
		answer[dna.T][4],
		answer[dna.T][5],
		answer[dna.T][6],
		answer[dna.T][7],
		answer[dna.T][8],
		answer[dna.T][9],
	)
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

// calculation probability of each substitution using row and column sums
func convertToProb(input [][]float64) [][]float64 {
	var currRowSum float64
	var currRow, currColumn int
	var output [][]float64 = make([][]float64, len(input))
	for currRow = range output {
		output[currRow] = make([]float64, len(input[currRow]))
	}
	for currRow = range input {
		currRowSum = 0
		for currColumn = range input[currRow] {
			currRowSum += input[currRow][currColumn]
		}
		for currColumn = range input[currRow] {
			output[currRow][currColumn] = input[currRow][currColumn] / currRowSum
		}
	}
	return output
}

// getPileCoverage calculates the coverage, defined as the number of As, Cs, Gs, and Ts observed at a position.
func getPileCoverage(p sam.Pile) int {
	var answer int
	answer += p.CountF[0] + p.CountR[0]
	answer += p.CountF[1] + p.CountR[1]
	answer += p.CountF[2] + p.CountR[2]
	answer += p.CountF[3] + p.CountR[3]
	return answer
}

// updateErrorEstimate updates the EmpiricalErrorEstimator for a given pileup. Epsilon is estimated as the empirical proportion
// of errors in homozygous A and T sites, whereas LambdaPlusEpsilon can be estimated from the empirical proportion of Ts in homozygous
// CC piles and As in homozygous GG piles.
func updateErrorEstimate(errorEstimator *EmpiricalErrorEstimator, baseCall sam.DiploidBase, p sam.Pile) {
	switch baseCall {
	case sam.AA:
		errorEstimator.NumEpsilon = errorEstimator.NumEpsilon + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
		errorEstimator.TotalEpsilon = errorEstimator.TotalEpsilon + p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	case sam.CC:
		errorEstimator.NumLambdaPlusEpsilon = errorEstimator.NumLambdaPlusEpsilon + p.CountF[dna.T] + p.CountR[dna.T]
		errorEstimator.TotalLambdaPlusEpsilon = errorEstimator.TotalLambdaPlusEpsilon + p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	case sam.GG:
		errorEstimator.NumLambdaPlusEpsilon = errorEstimator.NumLambdaPlusEpsilon + p.CountF[dna.A] + p.CountR[dna.A]
		errorEstimator.TotalLambdaPlusEpsilon = errorEstimator.TotalLambdaPlusEpsilon + p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	case sam.TT:
		errorEstimator.NumEpsilon = errorEstimator.NumEpsilon + p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G]
		errorEstimator.TotalEpsilon = errorEstimator.TotalEpsilon + p.CountF[dna.A] + p.CountF[dna.C] + p.CountF[dna.G] + p.CountF[dna.T] + p.CountR[dna.A] + p.CountR[dna.C] + p.CountR[dna.G] + p.CountR[dna.T]
	}
}
