package main

import (
	"flag"
	"fmt"
	"log"
	"os"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
)

// create struct for keeping list of arguments clean
type PriorSettings struct {
	SamFileName         string
	ReferenceFile       string
	OutFile             string
	Epsilon             float64
	LikelihoodCacheSize int
	PseudoCount         float64
	AsCounts            bool
}

// estimation of error rates
type EmpiricalErrorEstimator struct {
	NumEpsilon             int
	NumLambdaPlusEpsilon   int
	TotalEpsilon           int
	TotalLambdaPlusEpsilon int
}

// help message
func priorUsage(priorFlags *flag.FlagSet) {
	fmt.Print("samAssembler prior - Construct an empirical conditional Dirichlet prior for output diploid" +
		"genotypes based on maximum likelihood estimation of genotypes from aligned short reads.\n" +
		"Usage: \n" +
		"samAssembler [options] prior reads.sam/bam ref.fa output.txt\n" +
		"Options:\n")
	priorFlags.PrintDefaults()
}

// define and setup arguments for prior creation
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
	}
	SamAssemblerPrior(s)
}

// create prior for samAssembler pipeline
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
		// if the baseCall doesn't return a genotype with an unknown base, then add the base to the transition matrix
		if refBase < 4 {
			baseCall = sam.DiploidBaseCallFromPile(p, refBase, diploidBasePriorCache, homozygousCache, heterozygousCache, sam.AncientLikelihoodCache{}, s.Epsilon, 0)
			if baseCall < 10 {
				answer[refBase][baseCall]++
				updateErrorEstimate(&errorEstimator, baseCall, p)
			}
		}
	}

	EpsilonEstimate := float64(errorEstimator.NumEpsilon) / float64(errorEstimator.TotalEpsilon)
	LambdaEstimate := (float64(errorEstimator.NumLambdaPlusEpsilon) / float64(errorEstimator.TotalLambdaPlusEpsilon)) - EpsilonEstimate

	// convert matrix to probability
	if !s.AsCounts {
		answer = convertToProb(answer)
	}
	// format output file using aliases
	out := fileio.EasyCreate(s.OutFile)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Epsilon: %v\n", EpsilonEstimate)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Lambda: %v\n", LambdaEstimate)
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
