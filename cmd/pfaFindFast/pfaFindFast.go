// Command Group: "Sequence Evolution & Reconstruction"

// For two pFasta sequences, calculate a score at each base for the distance that separates the two pDNA bases, and return the average over a given window size
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
	InFile          string
	OutFile         string
	FirstQueryName  string
	SecondQueryName string
	WindowSize      int
	RefChromName    string
	RemoveN         bool
	LongOutput      bool
	DivergenceRate  float64
	OutputAlnPos    bool
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

	// Conversion to valid sequences. TODO: remove after debugging?
	firstQuery = pFasta.MakeValid(firstQuery)
	secondQuery = pFasta.MakeValid(secondQuery)

	// QC firstQuery and secondQuery sequences. TODO: remove after debugging?
	//pFasta.QC(firstQuery)
	//pFasta.QC(secondQuery)
	fmt.Printf("QC complete\n")
	//fmt.Printf("firstQuery: %v\n", firstQuery)
	//fmt.Printf("secondQuery: %v\n", secondQuery)
	speedyWindowDifference(reference, firstQuery, secondQuery, s)
}

// TODO: cmd for faToPfa?
// TODO: modify usage and flags as needed
func usage() {
	fmt.Print(
		"pfaFindFast - For two pFasta sequences, calculate a score at each base for the distance that separates the two pDNA bases, and return the average over a given window size\n" +
			"Usage:\n" +
			" pfaFindFast input.fa output.bed\n" +
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
	var longOutput *bool = flag.Bool("longOutput", false, "Print percent diverged and raw -Log10PValue in output. Requires the 'divergenceRate' argument.")
	var divergenceRate *float64 = flag.Float64("divergenceRate", math.MaxFloat64, "Set the null divergence rate for p value calculations with 'longOutput'.")
	var outputAlnPos *bool = flag.Bool("outputAlnPos", false, "Print the alignment position of the window's start in output as the last column.")
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
		InFile:          inFile,
		OutFile:         outFile,
		FirstQueryName:  *firstQueryName,
		SecondQueryName: *secondQueryName,
		WindowSize:      *windowSize,
		RefChromName:    *refChromName,
		RemoveN:         *removeN,
		LongOutput:      *longOutput,
		DivergenceRate:  *divergenceRate,
		OutputAlnPos:    *outputAlnPos,
	}

	pfaFindFast(s)
}
