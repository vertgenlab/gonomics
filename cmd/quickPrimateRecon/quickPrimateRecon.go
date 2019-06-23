package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func quickPrimateRecon(infile string, outfile string) {
	records := fasta.Read(infile)
	outputRecord := &fasta.Fasta{Name: "Human_Chimp_Ancestor"}

	var outputBase dna.Base	
	
	if len(records) != 5 {
		log.Fatalf("Wrong number of sequences, expecting five, found %d.\n", len(records))
	}

	//confirm alignment lengths are all the same
	firstLength := len(records[0].Seq)
	for i := 1; i < len(records); i++ {
		if len(records[i].Seq) != firstLength {
			log.Fatalf("Sequence %d is the wrong length.\n", i+1)
		}
	}

	var Ncount int = 0
	var allMatch int = 0
	var humanMatchChimporBonobo int = 0
	var humanChange int = 0
	var gorillaVote int = 0
	var orangutanVote int = 0
	var messyBase int = 0

	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]
	var humanLength int = len(human.Seq)

	for j := 0; j < len(human.Seq); j++ {
		outputBase = human.Seq[j]

		if !isThereABase(records, j) {
			outputBase = dna.Gap
		} else if human.Seq[j] == dna.N {
			Ncount++
			outputBase = human.Seq[j]
		} else if human.Seq[j] == chimp.Seq[j] && human.Seq[j] == bonobo.Seq[j] {
			outputBase = human.Seq[j]
			allMatch++
		} else if (human.Seq[j] == chimp.Seq[j] || human.Seq[j] == bonobo.Seq[j]) && human.Seq[j] != dna.Gap {
			outputBase = human.Seq[j]
			humanMatchChimporBonobo++
		} else if (chimp.Seq[j] == bonobo.Seq[j] && (chimp.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == orangutan.Seq[j])) && (chimp.Seq[j] != dna.N && chimp.Seq[j] != dna.Gap){
					outputBase = chimp.Seq[j]
					humanChange++
		} else if (human.Seq[j] == gorilla.Seq[j] || chimp.Seq[j] == gorilla.Seq[j] || bonobo.Seq[j] == gorilla.Seq[j]) && (gorilla.Seq[j] != dna.N && gorilla.Seq[j] != dna.Gap){
			outputBase = gorilla.Seq[j]
			gorillaVote++
		} else if (human.Seq[j] == orangutan.Seq[j] || chimp.Seq[j] == orangutan.Seq[j] || bonobo.Seq[j] == orangutan.Seq[j] || gorilla.Seq[j] == orangutan.Seq[j]) && (orangutan.Seq[j] != dna.N && orangutan.Seq[j] != dna.Gap) {
			outputBase = orangutan.Seq[j]
			orangutanVote++
		} else {
			if human.Seq[j] != dna.Gap {
				outputBase = human.Seq[j]
			} else {
				outputBase = dna.N
			}
			messyBase++
		}
		outputRecord.Seq = append(outputRecord.Seq, outputBase)
	}

	var outputSlice []*fasta.Fasta

	for k := 0; k < len(records); k++ {
		outputSlice = append(outputSlice, records[k])
	}

	outputSlice = append(outputSlice, outputRecord)

	fmt.Printf("Read %v bases, NCount: %v. allMatch: %v. humanMatchChimporBonobo: %v. humanChange: %v. gorillaVote: %v. orangutanVote: %v. messyBase: %v.\n",
			humanLength, Ncount, allMatch, humanMatchChimporBonobo, humanChange, gorillaVote, orangutanVote, messyBase)

	fasta.Write(outfile, outputSlice)
}

func isThereABase(records []*fasta.Fasta, j int) bool {
	human := records[0]
	bonobo := records[1]
	chimp := records[2]
	orangutan := records[3]
	gorilla := records[4]

	if (human.Seq[j] != dna.Gap || chimp.Seq[j] != dna.Gap || bonobo.Seq[j] != dna.Gap) && (gorilla.Seq[j] != dna.Gap || orangutan.Seq[j] != dna.Gap) {
			return true 
	} else if human.Seq[j] != dna.Gap && (chimp.Seq[j] != dna.Gap || bonobo.Seq[j] != dna.Gap) {
		return true
	}
	return false
}

func usage() {
	fmt.Print(
		"quickPrimateRecon - Returns maximum likelihood sequence from an HBCGO primate alignment\n" +
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
