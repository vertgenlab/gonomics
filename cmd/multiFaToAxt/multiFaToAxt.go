package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func multiFaToAxt(inFile string, rName string, qName string, outFile string) {
	var answer []axt.Axt = make([]axt.Axt, 0)
	var currAxt axt.Axt
	records := fasta.Read(inFile)
	if len(records) != 2 {
		log.Fatalf("multiFaToAxt accepts only a pairwise multiFa alignment. Expected 2 entries, found %d.\n", len(records))
	}
	//var recordsCopy []fasta.Fasta = fasta.RemoveGaps(records)

	currAxt = axt.Axt{
		RName:      rName,
		RStart:     0,
		//REnd:       len(recordsCopy[0].Seq),
		REnd:       len(records[0].Seq),
		QName:      qName,
		QStart:     0,
		//QEnd:       len(recordsCopy[1].Seq),
		QEnd:       len(records[1].Seq),
		QStrandPos: true,
		Score:      100, //dummy value, hardcoded
		RSeq:       records[0].Seq,
		QSeq:       records[1].Seq,
	}

	answer = append(answer, currAxt)
	axt.Write(outFile, answer)
}

func usage() {
	fmt.Print(
		"multiFaToAxt - Convert a multiFa format alignment to an AXT format alignment file.\n" +
			"Usage:\n" +
			" multiFaToAxt input.fa RName QName output.axt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 4
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	RName := flag.Arg(1)
	QName := flag.Arg(2)
	outFile := flag.Arg(3)

	multiFaToAxt(inFile, RName, QName, outFile)
}
