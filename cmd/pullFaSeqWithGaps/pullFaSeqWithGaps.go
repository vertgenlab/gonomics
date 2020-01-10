package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"strconv"
)

func pullFaSeqWithGaps (infile string, outfile string, start int64, end int64) {
	if !(start < end) {
		log.Fatalf("Invalid arguments, start must be lower than end")
	}
	records := fasta.Read(infile)
	var ans []*fasta.Fasta
	var refCounter int64 = 0
	var startCounter int64 = 0
	var endCounter int64 = 0

	for t := 0; refCounter < start; t++ {
		startCounter++
		if t == len(records[0].Seq) {
			log.Fatalf("Ran out of chromosome")
		} else if records[0].Seq[t] != dna.Gap {
			refCounter++
		}
	}

	refCounter = 0
	for n := 0; refCounter < end; n++ {
		endCounter++
		if n == len(records[0].Seq) {
			log.Fatalf("Ran off the chromosome")
		} else if records[0].Seq[n] != dna.Gap {
			refCounter++
		}
	}
	fmt.Printf("Start: %d. refCounter: %d. alignCounter: %d\n", start, refCounter, startCounter)

	for i := 0; i < len(records); i++ {
		ans = append(ans, fasta.Extract(records[i], startCounter, endCounter, records[i].Name))
	}
	fasta.Write(outfile, ans)
}

func usage() {
	fmt.Print(
		"multFaPullSeqs - Pull sub-sequence from multiple fasta alignment for each entry, retaining gaps and thus multiple alignment format.\n" +
		"Usage:\n" +
		"multFaPullSeqs mult.fa out.fa start end\n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	start, _ := strconv.ParseInt(flag.Arg(2), 10, 64)
	end, _ := strconv.ParseInt(flag.Arg(3), 10, 64)

	pullFaSeqWithGaps(infile, outfile, start, end)
}