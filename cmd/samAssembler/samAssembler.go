// Command Group: "Data Conversion"

// Reference-based diploid assembly of aligned short reads
package main

import (
	"flag"
	"fmt"
	"log"
)

const bufferSize = 10_000_000

func usage() {
	fmt.Print(
		"samAssembler - Reference-based diploid assembly of aligned short reads.\n" +
			"Can be used in three modes:\n\t'build' generates diploid assemblies.\n" +
			"\t'score' validates assembly accuracy with a five way alignment including the known divergent sequences\n" +
			"\t'prior' constructs an empirical prior for output diploid genotypes based on a maximum likelihood estimate from the input reads.\n" +
			"Enter: 'samAssembler build' OR 'samAssembler score' OR 'samAssembler prior' to view usage and options.\n")
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Expecting at least one argument, received 0.")
	}

	var mode = flag.Arg(0)

	switch mode {
	case "build":
		var delta *float64 = flag.Float64("delta", 0.01, "Set the expected divergence frequency.")
		var gamma *float64 = flag.Float64("gamma", 3, "Set the expected transition bias.")
		var epsilon *float64 = flag.Float64("epsilon", 0.01, "Set the expected misclassification error rate.")
		var kappa *float64 = flag.Float64("kappa", 0.1, "Set the expected proportion of divergent sites that are INDELs.")
		var lambda *float64 = flag.Float64("lambda", 0, "Set the expected rate of cytosine deamination.")
		var multiFaDir *string = flag.String("multiFaDir", "", "Output the reference and generated sequences as an aligned multiFa, each file by chrom.")
		var qNameA *string = flag.String("qNameA", "QueryA", "Set the qName for the first generated chromosome in the optional multiFa output.")
		var qNameB *string = flag.String("qNameB", "QueryB", "Set the qName for the second generated chromosome in the optional multiFa output.")
		var likelihoodCacheSize *int = flag.Int("likelihoodCacheSize", 100, "Set the maximum dimension of the likelihood caches. Should be slightly larger than highest expected pile depth.")
		var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
		var verbose *int = flag.Int("verbose", 0, "Set to 1 to enable additional debug prints.")
		var flatPrior *bool = flag.Bool("flatPrior", false, "Use a flat prior instead of the default informative prior distribution.")
		var empricalPrior *string = flag.String("empiricalPrior", "", "Use an empirical prior instead, based on an input prior file. New empirical priors can be generated with the 'prior' subcommand.")
		var expectedNumArgs int = 5

		flag.Usage = buildUsage
		log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
		flag.Parse()

		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}

		inFile := flag.Arg(1)
		refFile := flag.Arg(2)
		outFileA := flag.Arg(3)
		outFileB := flag.Arg(4)

		s := BuildSettings{
			SamFileName:         inFile,
			RefFile:             refFile,
			OutFileA:            outFileA,
			OutFileB:            outFileB,
			MultiFaDir:          *multiFaDir,
			qNameA:              *qNameA,
			qNameB:              *qNameB,
			Delta:               *delta,
			Gamma:               *gamma,
			Epsilon:             *epsilon,
			Kappa:               *kappa,
			Lambda:              *lambda,
			LikelihoodCacheSize: *likelihoodCacheSize,
			SetSeed:             *setSeed,
			Verbose:             *verbose,
			FlatPrior:           *flatPrior,
			EmpiricalPrior:      *empricalPrior,
		}

		samAssemblerBuild(s)
	case "score":
		var expectedNumArgs int = 4
		flag.Usage = scoreUsage
		log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
		flag.Parse()

		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		scoreType := flag.Arg(1)
		inFileList := flag.Arg(2)
		outFile := flag.Arg(3)
		s := ScoreSettings{
			ScoreType:  scoreType,
			InFileList: inFileList,
			OutFile:    outFile,
		}
		samAssemblerScore(s)
	case "prior":
		var epsilon *float64 = flag.Float64("epsilon", 0.01, "Set the expected misclassification error rate.")
		var likelihoodCacheSize *int = flag.Int("likelihoodCacheSize", 100, "Set the maximum dimension of the likelihood caches. Should be slightly larger than highest expected pile depth.")
		var pseudoCount *float64 = flag.Float64("pseudoCount", 0.01, "Prime the empirical prior with pseudoCounts to avoid prior probabilities of 0.\n"+
			"Increasing this value regularizes the prior against overfitting.\n"+
			"The model will underfit, biased towards Jukes-Cantor, if this value is set too high.")
		var asCounts *bool = flag.Bool("asCounts", false, "Return output as counts, instead of probabilities. Useful for debugging, but count matrices cannot be used directly in 'build'.")
		var expectedNumArgs int = 4
		flag.Usage = priorUsage
		log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
		flag.Parse()

		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}

		samFileName := flag.Arg(1)
		referenceFileName := flag.Arg(2)
		outFileName := flag.Arg(3)

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
	default:
		log.Fatalf("Unknown mode. samAssembler can be run with the first argument as 'build' or 'score'.")
	}
}
