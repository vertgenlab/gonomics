package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faFormat(inFile string, outFile string, lineLength int) {
	records := fasta.Read(inFile)

	file := fileio.MustCreate(outFile)
	defer file.Close()

	fasta.WriteToFileHandle(file, records, lineLength)
}

func usage() {
	fmt.Print(
		"faFormat - reformat the sequences in a fasta file\n" +
			"Usage:\n" +
			" faFormat input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var lineLength *int = flag.Int("lineLength", 50, "wrap sequence lines after this many characters")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFormat(inFile, outFile, *lineLength)
}
