package main

import (
	"fmt"
	"log"
	"flag"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/dna"
)

func faBaseCounts(infile string) {
	records := fasta.Read(infile)
	A := 0
	T := 0
	C := 0
	G := 0
	N := 0
	Gap := 0

	for i := 0; i < len(records); i++ {
		for k := 0; k < len(records[i].Seq); k++ {
			if records[i].Seq[k] == dna.A {
				A++
			} else if records[i].Seq[k] == dna.C {
				C++
			} else if records[i].Seq[k] == dna.T {
				T++
			} else if records[i].Seq[k] == dna.G {
				G++
			} else if records[i].Seq[k] == dna.Gap {
				Gap++
			} else if records[i].Seq[k] == dna.N {
				N++
			}
		}
		fmt.Printf("Name: %s. A: %d. T: %d. C: %d. G: %d. Gap: %d. N: %d.\n", records[i].Name, A, T, C, G, Gap, N)
	}
}

func usage() {
	fmt.Print(
		"faBaseCounts - Returns the counts for each base.\n" +
		"Usage:\n" +
		" faBaseCounts infile.fa\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)

	faBaseCounts(inFile)
}