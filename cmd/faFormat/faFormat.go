package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

type Settings struct {
	InFile string
	OutFile string
	LineLength int
	TrimName bool
	ToUpper bool
	RevComp bool
	NoGaps bool
}

func faFormat(s Settings) {
	records := fasta.Read(s.InFile)

	if s.NoGaps {
		fasta.RemoveGaps(records)
	}

	for i := range records {
		if s.TrimName {
			records[i] = fasta.TrimName(records[i])
		}
		if s.ToUpper {
			fasta.ToUpper(records[i])
		}
		if s.RevComp {
			fasta.ReverseComplement(records[i])
		}
	}

	file := fileio.EasyCreate(s.OutFile)
	defer file.Close()

	fasta.WriteToFileHandle(file, records, s.LineLength)
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
	var trimName *bool = flag.Bool("trimName", false, "if a fasta name contains spaces, retains only the first space delimited field")
	var toUpper *bool = flag.Bool("toUpper", false, "Convert all DNA bases to upper case.")
	var revComp *bool = flag.Bool("revComp", false, "Return the reverse complement for each sequence.")
	var noGaps *bool = flag.Bool("noGaps", false, "Remove gaps from all input sequences.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		InFile: inFile,
		OutFile: outFile,
		LineLength: *lineLength,
		TrimName: *trimName,
		RevComp: *revComp,
		ToUpper: *toUpper,
		NoGaps: *noGaps,
	}

	faFormat(s)
}
