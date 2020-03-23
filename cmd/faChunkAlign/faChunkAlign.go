package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"os"
	"strconv"
)

func faChunkAlign(inFile string, chunkSize int, gapOpen int64, gapExtend int64, outFile string) {
	log.Printf("Reading %s...\n", inFile)
	records := fasta.Read(inFile)
	log.Printf("Successfully read %d sequences from fasta file.\n", len(records))

	log.Printf("Aligning sequences...")
	records = align.AllSeqAffineChunk(records, align.HumanChimpTwoScoreMatrix, gapOpen, gapExtend, chunkSize)

	log.Printf("Writing aligned sequences to %s...", outFile)
	fasta.Write(outFile, records)
	log.Print("Done\n")
}

func usage() {
	fmt.Fprint(os.Stderr,
		"faChunkAlign - Align two or more sequences by \"chunks\" of bases\n"+
			"                instead of by single bases.  Each sequence must\n"+
			"                have a length that is divisible by the chunk size.\n"+
			"Usage:\n"+
			" faChunkAlign unaligned.fa chunkSize aligned.fa\n"+
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var gapOpen *int64 = flag.Int64("gapOpen", 300, "Penalty for opening a gap")
	var gapExtend *int64 = flag.Int64("gapExtend", 40, "Penalty for extending a gap")
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	chunkSize, err := strconv.Atoi(flag.Arg(1))
	common.ExitIfError(err)
	outFile := flag.Arg(2)

	faChunkAlign(inFile, chunkSize, *gapOpen*-1, *gapExtend*-1, outFile)
}
