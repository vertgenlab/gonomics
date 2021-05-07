package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"log"
)

func randSeq(outFile string, GC float64, numSeq int, lenSeq int, randSeed bool, setSeed int64) {
	common.RngSeed(randSeed, setSeed)
	file := fileio.EasyCreate(outFile)
	defer file.Close()
	for i := 0; i < numSeq; i++ {
		fasta.WriteFasta(file, fasta.Fasta{Name: fmt.Sprintf("Sequence_%v", i), Seq: simulate.RandIntergenicSeq(GC, lenSeq)}, 50)
	}
}

func usage() {
	fmt.Print(
		"randSeq - returns pseudorandomly generated DNA sequences in fasta format.\n" +
			"Usage:\n" +
			" randSeq out.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var GC *float64 = flag.Float64("GC", 0.41, "Set the expected value of the GC content of the output sequences.")
	var numSeq *int = flag.Int("numSeq", 10, "Set the number of sequences to generate.")
	var lenSeq *int = flag.Int("lenSeq", 500, "Set the length of sequences to generate.")
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	outFile := flag.Arg(0)
	randSeq(outFile, *GC, *numSeq, *lenSeq, *randSeed, *setSeed)
}