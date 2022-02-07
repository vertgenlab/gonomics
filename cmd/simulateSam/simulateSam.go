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
)

func usage() {
	fmt.Print(
		"simulateSam - Simulate alignments to a reference sequence\n\n" +
			"Currently only generates illumina-style paired sequencing data" +
			"Read Length : 150bp, Average Insert Size: 50bp" +
			"Usage:\n" +
			" simulateSam [options] -i ref.fasta\n" +
			"options:\n")
	flag.PrintDefaults()
}

func simulateSam(refFile, outFile string, numReads int, bamOutput bool) {
	ref := fasta.Read(refFile)
	out := fileio.EasyCreate(outFile)
	header := sam.GenerateHeader(fasta.ToChromInfo(ref), nil, sam.Unsorted, sam.None)

	var bw *sam.BamWriter
	if bamOutput {
		bw = sam.NewBamWriter(out, header)
	} else {
		sam.WriteHeaderToFileHandle(out, header)
	}

	var readsPerContig []int = getReadsPerContig(ref, numReads)
	var reads []sam.Sam
	for i := range ref {
		reads = simulate.IlluminaSam(ref[i].Name, ref[i].Seq, readsPerContig[i])
		if bamOutput {
			for i := range reads {
				sam.WriteToBamFileHandle(bw, reads[i], 0)
			}
		} else {
			for i := range reads {
				sam.WriteToFileHandle(out, reads[i])
			}
		}
	}

	var err error
	if bamOutput {
		err = bw.Close()
		exception.PanicOnErr(err)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func getReadsPerContig(ref []fasta.Fasta, numReads int) []int {
	// get probability weights for each contig based on length
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

func main() {
	outFile := flag.String("o", "stdout", "output file")
	refFile := flag.String("r", "", "reference fasta file")
	numReads := flag.Int("n", 100, "number of read pairs to generate")
	bamOutput := flag.Bool("b", false, "output file as .bam instead of .sam")
	seed := flag.Int64("seed", 1, "set the seed for the simulation")
	flag.Usage = usage
	flag.Parse()

	if *refFile == "" {
		flag.Usage()
		log.Fatal("ERROR: must input a reference fasta file (-r)")
	}

	rand.Seed(*seed)
	simulateSam(*refFile, *outFile, *numReads, *bamOutput)
}
