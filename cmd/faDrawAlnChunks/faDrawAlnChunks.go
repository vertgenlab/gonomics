// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"image/png"
	"log"
	"os"
	"strconv"
)

func faDrawAlnChunks(faAlnFilename string, chunkSize int, imageOutFilename string) {
	aln := fasta.Read(faAlnFilename)
	img, err := align.DrawAlignedChunks(aln, chunkSize, 6, 12)
	exception.PanicOnErr(err)
	imgOutFile, err := os.Create(imageOutFilename)
	exception.PanicOnErr(err)
	err = png.Encode(imgOutFile, img)
	exception.PanicOnErr(err)
	err = imgOutFile.Close()
	exception.PanicOnErr(err)
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
	exception.PanicOnErr(err)
	imageOutFilename := flag.Arg(2)

	faDrawAlnChunks(alignmentFilename, chunkSize, imageOutFilename)
}
