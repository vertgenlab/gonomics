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

	var totalLen int
	for i := range ref {
		totalLen += len(ref[i].Seq)
	}

	var proportionOfReads float64
	var readsPerContig int
	var reads []sam.Sam
	for i := range ref {
		proportionOfReads = float64(len(ref[i].Seq)) / float64(totalLen)
		readsPerContig = int(float64(numReads) * proportionOfReads)
		reads = simulate.IlluminaSam(ref[i].Name, ref[i].Seq, readsPerContig)

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
