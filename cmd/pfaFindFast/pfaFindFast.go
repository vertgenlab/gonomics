// Command Group: "Sequence Evolution & Reconstruction"

// For two pFasta sequences, calculate a score at each base for the divergence between the two pDNA bases, convert the score to a binary substitution/non-substitution call, and return the calculations over a given window size
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"log"
	"math"
	"os"
	"runtime/pprof"
)

type Settings struct {
	InFile                  string
	OutFile                 string
	FirstQueryName          string
	SecondQueryName         string
	WindowSize              int
	RefChromName            string
	RemoveN                 bool
	LongOutput              bool
	DivergenceRate          float64
	OutputAlnPos            bool
	BaseDistToDivThreshold  float64
	BaseDotToSubstThreshold float64
	ConfidentThreshold      float32
}

func pfaFindFast(s Settings) {
	var reference, firstQuery, secondQuery []pDna.Float32Base
	var found bool

	records := pFasta.Read(s.InFile)    // fasta.Read already makes sure that sequence names are unique
	recordsMap := pFasta.ToMap(records) // generate map[string][]dna.Base. No order but can do key-value search

	if len(records) < 2 {
		log.Fatalf("Error: There must be at least 2 pFasta records in the input file.\n")
	}

	//if reference and query names were not specified in the command, then use 1st and 2nd fasta records' names, respectively
	//allows backward compatibility with previously-written code
	if s.FirstQueryName == "" {
		firstQuery = records[0].Seq
	} else {
		_, found = recordsMap[s.FirstQueryName]
		if found {
			firstQuery = recordsMap[s.FirstQueryName]
		} else {
			log.Fatalf("Error: first query name is not found in the input file.\n")
		}
	}
	if s.SecondQueryName == "" {
		secondQuery = records[1].Seq
	} else {
		_, found = recordsMap[s.SecondQueryName]
		if found {
			secondQuery = recordsMap[s.SecondQueryName]
		} else {
			log.Fatalf("Error: second query name is not found in the input file.\n")
		}
	}

	reference = records[0].Seq // reference will always be the 1st pFasta record in a multi-pFa input file

	if !(len(reference) == len(firstQuery) && len(reference) == len(secondQuery)) {
		log.Fatalf("Error: Reference, first query, and second query sequences are not all of equal length.\n")
	}

	speedyWindowDifference(reference, firstQuery, secondQuery, s)
}

func usage() {
	fmt.Print(
		"pfaFindFast - For two pFasta sequences, calculate scores at each base for the divergence that separates the two pDNA bases, and return a bed file with calculations.\n" +
			"Usage:\n" +
			" pfaFindFast input.pFa output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var firstQueryName *string = flag.String("firstQueryName", "", "Specify the name of the first query sequence")
	var secondQueryName *string = flag.String("secondQueryName", "", "Specify the name of the second query sequence")
	var windowSize *int = flag.Int("windowSize", 1000, "Specify the window size")
	var refChromName *string = flag.String("chrom", "", "Specify a chrom name of the reference sequence")
	var removeN *bool = flag.Bool("removeN", false, "Excludes bed regions with Ns in the reference from the output.")
	var longOutput *bool = flag.Bool("longOutput", false, "Print these additional fields in output: percent diverged, raw -Log10PValue, total number of substitutions in the window, total number of gaps in the window, window size in the reference sequence, total number of confident bases for the second query sequence in the window. Requires the 'divergenceRate' argument.")
	var divergenceRate *float64 = flag.Float64("divergenceRate", math.MaxFloat64, "Set the null divergence rate for p value calculations with 'longOutput'.")
	var outputAlnPos *bool = flag.Bool("outputAlnPos", false, "Print this additional field in output: the alignment position of the window's start.")
	var baseDotToSubstThreshold *float64 = flag.Float64("baseDotToSubstThreshold", 0.8, "Specify the threshold above which to call the 1 - dot product of 2 bases a substitution. Defaults to 0.8.")
	var baseDistToDivThreshold *float64 = flag.Float64("baseDistToDivThreshold", 0.7, "Specify the threshold above which to call the distance between 2 bases a divergence. Defaults to 0.7.")
	var confidentThreshold *float64 = flag.Float64("confidentThreshold", 0.8, "Specify the threshold above which to call the most likely base at a position confident. Defaults to 0.8.")
	// for go proflier cpu
	var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

	if *longOutput && *divergenceRate == math.MaxFloat64 {
		log.Fatalf("Error: must set a 'divergenceRate' if using the 'longOutput' option.\n")
	}

	if *divergenceRate != math.MaxFloat64 {
		if *divergenceRate < 0 || *divergenceRate > 1 {
			log.Fatalf("Error: divergence rate must be a value between 0 and 1.\n")
		}
	}

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	// for go profiler cpu
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d.\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		InFile:                  inFile,
		OutFile:                 outFile,
		FirstQueryName:          *firstQueryName,
		SecondQueryName:         *secondQueryName,
		WindowSize:              *windowSize,
		RefChromName:            *refChromName,
		RemoveN:                 *removeN,
		LongOutput:              *longOutput,
		DivergenceRate:          *divergenceRate,
		OutputAlnPos:            *outputAlnPos,
		BaseDistToDivThreshold:  *baseDistToDivThreshold,
		BaseDotToSubstThreshold: *baseDotToSubstThreshold,
		ConfidentThreshold:      float32(*confidentThreshold),
	}

	pfaFindFast(s)
}
