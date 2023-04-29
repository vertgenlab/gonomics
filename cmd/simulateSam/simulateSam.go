// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math/rand"
	"sort"
	"strings"
)

type Settings struct {
	OutFile                 string
	RefFile                 string
	NumReads                int
	Coverage                float64
	ReadLength              int
	FlatError               float64
	FragmentLength          int
	FragmentStdDev          float64
	AncientErrorRate        float64
	GeometricParam          float64
	SetSeed                 int64
	DeaminationDistribution string
}

func simulateSam(s Settings) {
	rand.Seed(s.SetSeed)
	ref := fasta.Read(s.RefFile)
	out := fileio.EasyCreate(s.OutFile)
	header := sam.GenerateHeader(fasta.ToChromInfo(ref), nil, sam.Unsorted, sam.None)

	var bw *sam.BamWriter
	var bamOutput bool
	if strings.HasSuffix(s.OutFile, ".bam") {
		bw = sam.NewBamWriter(out, header)
		bamOutput = true
	} else {
		sam.WriteHeaderToFileHandle(out, header)
		bamOutput = false
	}

	flatBinomialAlias := numbers.MakeBinomialAlias(s.ReadLength, s.FlatError)
	ancientBinomialAlias := numbers.MakeBinomialAlias(s.ReadLength, s.AncientErrorRate)

	var deaminationDistributionSlice []int = make([]int, s.ReadLength)

	var readsPerContig = getReadsPerContig(ref, s.NumReads, s.Coverage, s.ReadLength)
	for i := range ref {
		simulate.IlluminaPairedSam(ref[i].Name, ref[i].Seq, readsPerContig[i], s.ReadLength, s.FragmentLength, s.FragmentStdDev, s.FlatError, s.AncientErrorRate, flatBinomialAlias, ancientBinomialAlias, s.GeometricParam, out, bw, bamOutput, deaminationDistributionSlice)
	}

	var err error

	if s.DeaminationDistribution != "" {
		deaminationOut := fileio.EasyCreate(s.DeaminationDistribution)
		_, err = fmt.Fprintf(deaminationOut, "Position\tCount\n")
		for i := range deaminationDistributionSlice {
			_, err = fmt.Fprintf(deaminationOut, "%v\t%v\n", i, deaminationDistributionSlice[i])
			exception.PanicOnErr(err)
		}
		err = deaminationOut.Close()
		exception.PanicOnErr(err)
	}

	if bamOutput {
		err = bw.Close()
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// get probability weights for each contig based on length.
func getReadsPerContig(ref []fasta.Fasta, numReads int, coverage float64, readLen int) []int {
	var totalLen int
	for i := range ref {
		totalLen += len(ref[i].Seq)
	}
	contigWeights := make([]float64, len(ref))
	for i := range ref {
		contigWeights[i] = float64(len(ref[i].Seq)) / float64(totalLen)
	}

	// make cumulative distribution
	cdf := make([]float64, len(contigWeights))
	cdf[0] = contigWeights[0]
	for i := 1; i < len(contigWeights); i++ {
		cdf[i] = cdf[i-1] + contigWeights[i]
	}

	//if the user specified coverage, we will set numReads.
	// readLen*2 for paired end gets bases per molecule.
	if coverage > 0 {
		numReads = int(coverage * float64(totalLen) / float64(readLen*2))
	}

	// make random pulls from cumulative distribution
	readsPerContig := make([]int, len(ref))
	var chosenContig int
	var val float64
	for i := 0; i < numReads; i++ {
		val = rand.Float64()
		// binary search for smallest index with cumulative sum > random value
		chosenContig = sort.Search(len(cdf), func(i int) bool { return cdf[i] > val })
		readsPerContig[chosenContig]++
	}
	return readsPerContig
}

func usage() {
	fmt.Print(
		"simulateSam - Simulate alignments to a reference sequence\n\n" +
			"Generates illumina-style paired sequencing data." +
			"Usage:\n" +
			"simulateSam [options] ref.fasta out.sam/bam\n" +
			"options:\n")
	flag.PrintDefaults()
}

func preCheck(s Settings) {
	if s.Coverage < 0 {
		log.Fatalf("Error: coverage must be a value greater than 0. Found: %v.\n", s.Coverage)
	}
	if s.NumReads < 1 && s.Coverage <= 0 {
		log.Fatalf("Error: number of reads must be greater than 0. Found: %v.\n", s.NumReads)
	}
	if s.ReadLength < 1 {
		log.Fatalf("Error: read length must be greater than 0. Found: %v.\n", s.ReadLength)
	}
	if s.FragmentLength < 1 {
		log.Fatalf("Error: fragment length must be greater than 0. Found: %v.\n", s.FragmentLength)
	}
	if s.FragmentStdDev <= 0 {
		log.Fatalf("Error: fragment standard deviation must be greater than 0. Found: %v.\n", s.FragmentStdDev)
	}
	if s.FlatError < 0 || s.FlatError > 1 {
		log.Fatalf("Error: flat error rate must be a value between 0 and 1. Found: %v.\n", s.FlatError)
	}
	if s.GeometricParam < 0 || s.GeometricParam > 1 {
		log.Fatalf("Error: geometric distribution parameter must be a value between 0 and 1. Found: %v.\n", s.GeometricParam)
	}
	if s.AncientErrorRate < 0 || s.AncientErrorRate > 1 {
		log.Fatalf("Error: ancient error rate must be a value between 0 and 1. Found: %v.\n", s.AncientErrorRate)
	}
	if s.Coverage > 0 && s.NumReads > 0 {
		log.Fatalf("Error: User must specify either -coverage or -n, not both.")
	}
}

func main() {
	var expectedNumArgs int = 2
	numReads := flag.Int("n", 0, "number of read pairs to generate")
	coverage := flag.Float64("coverage", 0, "Set an expected coverage value instead of specifying the number of reads.")
	setSeed := flag.Int64("setSeed", 1, "set the seed for the simulation")
	readLength := flag.Int("readLength", 150, "Set the read length for each paired end.")
	fragmentLength := flag.Int("fragmentLength", 400, "Set the average library fragment size. Note that the minimum fragment length will be forced to be equal to readLength.")
	fragmentStdDev := flag.Float64("fragmentStdDev", 50, "Set the library fragment size standard deviation.")
	flatError := flag.Float64("flatErrorRate", 0, "Sets an error rate for bases in synthetic reads. Each base will appear as one of the three other bases in the generated read with this probability.")
	geometricParam := flag.Float64("geometricParam", 0.25, "For ancient DNA libraries, set the success probability parameter for the geometric distribution defining the spatial distribution of cytosine deamination.")
	ancientErrorRate := flag.Float64("ancientErrorRate", 0, "Sets the expected rate of cytosine deamination. p(C->T|C or G->A|G).")
	deaminationDistribution := flag.String("deaminationDistribution", "", "Specify a file in which to print the spatial distribution of cytosine deamination events.")
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	refFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		RefFile:                 refFile,
		OutFile:                 outFile,
		NumReads:                *numReads,
		Coverage:                *coverage,
		ReadLength:              *readLength,
		FlatError:               *flatError,
		FragmentLength:          *fragmentLength,
		FragmentStdDev:          *fragmentStdDev,
		GeometricParam:          *geometricParam,
		AncientErrorRate:        *ancientErrorRate,
		SetSeed:                 *setSeed,
		DeaminationDistribution: *deaminationDistribution,
	}

	preCheck(s)

	simulateSam(s)
}
