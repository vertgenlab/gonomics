// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
	"math/rand"
	"sort"
	"strings"
)

type Settings struct {
	OutFile        string
	RefFile        string
	NumReads       int
	ReadLength     int
	FragmentLength int
	FragmentStdDev float64
	SetSeed        int64
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

	var readsPerContig = getReadsPerContig(ref, s.NumReads)
	for i := range ref {
		simulate.IlluminaPairedSam(ref[i].Name, ref[i].Seq, readsPerContig[i], s.ReadLength, s.FragmentLength, s.FragmentStdDev, out, bw, bamOutput)
	}

	var err error
	if bamOutput {
		err = bw.Close()
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// get probability weights for each contig based on length
func getReadsPerContig(ref []fasta.Fasta, numReads int) []int {
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
			"Currently only generates illumina-style paired sequencing data" +
			"Read Length : 150bp, Average Insert Size: 50bp" +
			"Usage:\n" +
			" simulateSam [options] ref.fasta out.sam/bam\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	numReads := flag.Int("n", 100, "number of read pairs to generate")
	setSeed := flag.Int64("setSeed", 1, "set the seed for the simulation")
	readLength := flag.Int("readLength", 150, "Set the read length for each paired end.")
	fragmentLength := flag.Int("fragmentLength", 400, "Set the average library fragment size.")
	fragmentStdDev := flag.Float64("fragmentStdDev", 50, "Set the library fragment size standard deviation.")
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
		RefFile:        refFile,
		OutFile:        outFile,
		NumReads:       *numReads,
		ReadLength:     *readLength,
		FragmentLength: *fragmentLength,
		FragmentStdDev: *fragmentStdDev,
		SetSeed:        *setSeed,
	}

	simulateSam(s)
}
