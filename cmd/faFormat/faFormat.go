package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faFormat(inFile string, outFile string, lineLength int, trimName bool) {
	records := fasta.Read(inFile)

	for i := range records{
		if trimName {
			records[i] = fasta.TrimName(records[i])
		}

	}

	file := fileio.EasyCreate(outFile)
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
	var trimName *bool = flag.Bool("trimName", false, "if a fasta name contains spaces, retains only the first space delimited field" )
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFormat(inFile, outFile, *lineLength, *trimName)
}
