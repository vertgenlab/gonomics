package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func quickPrimateRecon(infile string, outfile string) {
	records, err := fasta.Read(infile)
	var outputRecord fasta.Fasta
	outputRecord.Name = "Human_Chimp_Ancestor"
	var outputBase dna.Base

	if err != nil {
		log.Fatal(err)
	}

	if len(records) != 6 {
		log.Fatalf("Wrong number of sequences, expecting six, found %d.\n", len(records))
	}

	//confirm alignment lengths are all the same
	firstLength := len(records[0].Seq)
	for i := 1; i < len(records); i++ {
		if len(records[i].Seq) != firstLength {
			log.Fatalf("Sequence %d is the wrong length.\n", i+1)
		}
	}

	human := &records[0]
	bonobo := &records[1]
	chimp := &records[2]
	gorilla := &records[3]
	orangutan := &records[4]
	rhesus := &records[5]

	for j := 0; j < len(human.Seq); j++ {
		outputBase = human.Seq[j]
		if human.Seq[j] != chimp.Seq[j] && chimp.Seq[j] == bonobo.Seq[j] {
			if chimp.Seq[j] == gorilla.Seq[j] {
				outputBase = chimp.Seq[j]
			} else if chimp.Seq[j] == orangutan.Seq[j] {
				outputBase = chimp.Seq[j]
			} else if chimp.Seq[j] == rhesus.Seq[j] {
				outputBase = chimp.Seq[j]
			}
		}
		outputRecord.Seq = append(outputRecord.Seq, outputBase)
	}

	var outputSlice []fasta.Fasta

	for k := 0; k < len(records); k++ {
		outputSlice = append(outputSlice, records[k])
	}

	outputSlice = append(outputSlice, outputRecord)

	fasta.Write(outfile, outputSlice)
}

func usage() {
	fmt.Print(
		"quickPrimateRecon - Returns maximum likelihood sequence from an HBC GOQ primate alignment\n" +
			"Usage:\n" +
			"  quickPrimateRecon input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	quickPrimateRecon(inFile, outFile)
}
