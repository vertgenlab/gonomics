package main

import (
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"time"
)

func usage() {
	fmt.Print(
		"mkDictionaryFa - indexing fasta with seeds to align to a reference \n
		Takes one input fasta file")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}


	ref := fasta.Read(flag.Arg(0))
	start := time.Now()
	ham5 := qDna.IndexRefSlidingWindow(qDna.FromFastaSlice(ref), 20)
	qDna.Write("gasAcu.ham5", ham5)
	fmt.Println("Time to index reference and write to ham5: ", t1.Sub(start))
}
