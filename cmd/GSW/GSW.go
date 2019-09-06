package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/qDna"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"runtime"
	"time"
)

func usage() {
	fmt.Print(
		"GSW aligner \nTakes 4 inputs:\n1) Reference Fasta 2) fastq 3) index 4) sam output name\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	ref := fasta.Read(flag.Arg(0))
	fq := fastq.Read(flag.Arg(1))
	start := time.Now()

	//Read in index file

	samAlignment := qDna.GSW(qDna.FromFastaSlice(ref), fq, qDna.Read(flag.Arg(2)))
	t1 := time.Now()
	header := sam.AlignmentHeader(ref)
	fmt.Println("Number of reads input: ", len(fq))
	fmt.Println("CPUs: ", runtime.NumCPU())
	fmt.Println("Time to align fastq file to reference: ", t1.Sub(start))
	samFile := sam.Sam{Header: header, Aln: samAlignment}
	fileName := flag.Arg(3) + ".sam"
	err := sam.Write(fileName, &samFile)
	if err != nil {
		log.Fatalf("Reading %s gave an error")
	}

}
