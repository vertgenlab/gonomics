package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"image/png"
	"log"
	"os"
	"strconv"
)

func faDrawAlnChunks(faAlnFilename string, chunkSize int, imageOutFilename string) {
	log.Printf("Reading %s...\n", faAlnFilename)
	aln, err := fasta.Read(faAlnFilename)
	common.ExitIfError(err)
	log.Printf("Successfully read %d sequences from fasta file.\n", len(aln))

	img, err := align.DrawAlignedChunks(aln, chunkSize, 6, 12)
	common.ExitIfError(err)

	log.Printf("Writing aligned sequences to %s...", imageOutFilename)
	imgOutFile, err := os.Create(imageOutFilename)
	common.ExitIfError(err)
	defer imgOutFile.Close()
	err = png.Encode(imgOutFile, img)
	common.ExitIfError(err)
	log.Print("Done\n")
}

func usage() {
	fmt.Fprint(os.Stderr,
		"faDrawAlnChunks - Align two or more sequeces by \"chunks\" of bases\n"+
			"                   instead of by single bases.  Each sequence must\n"+
			"                   have a length that is divisible by the chunk size.\n"+
			"Usage:\n"+
			" faDrawAlnChunks aligned.fa chunkSize imageOut.png\n"+
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	alignmentFilename := flag.Arg(0)
	chunkSize, err := strconv.Atoi(flag.Arg(1))
	common.ExitIfError(err)
	imageOutFilename := flag.Arg(2)

	faDrawAlnChunks(alignmentFilename, chunkSize, imageOutFilename)
}
