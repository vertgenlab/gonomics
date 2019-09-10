package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/qDna"
	"log"
	"time"
)

func usage() {
	fmt.Print(
		"mkDictionaryFa - indexing fasta with seeds to align to a reference \nTakes two inputs: 1) fasta file 2) Name of index\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
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
	qDna.Write(flag.Arg(1), ham5)
	t1 := time.Now()
	fmt.Println("Time to index reference and write to ham5: ", t1.Sub(start))
}
