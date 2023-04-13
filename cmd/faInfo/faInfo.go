// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func faInfo(inFile string, outFile string) {
	records := fasta.Read(inFile)
	out := fileio.EasyCreate(outFile)
	recordA, recordT, recordC, recordG, recordGap, recordN := 0, 0, 0, 0, 0, 0
	totalA, totalT, totalC, totalG, totalGap, totalN := 0, 0, 0, 0, 0, 0
	var err error

	_, err = fmt.Fprintf(out, "#Name\tA\tT\tC\tG\tGap\tN\n") // create the header line.
	exception.PanicOnErr(err)

	for i := range records {
		recordA, recordT, recordC, recordG, recordGap, recordN = 0, 0, 0, 0, 0, 0
		for k := range records[i].Seq {
			if records[i].Seq[k] == dna.A {
				recordA++
				totalA++
			} else if records[i].Seq[k] == dna.T {
				recordT++
				totalT++
			} else if records[i].Seq[k] == dna.C {
				recordC++
				totalC++
			} else if records[i].Seq[k] == dna.G {
				recordG++
				totalG++
			} else if records[i].Seq[k] == dna.Gap {
				recordGap++
				totalGap++
			} else if records[i].Seq[k] == dna.N {
				recordN++
				totalN++
			} else {
				log.Fatalf("Character '%s' encountered in the fasta. This is an illegal character.", dna.BaseToString(records[i].Seq[k]))
			}
		}
		_, err = fmt.Fprintf(out, "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", records[i].Name, recordA, recordT, recordC, recordG, recordGap, recordN)
		exception.PanicOnErr(err)
	}
	_, err = fmt.Fprintf(out, "total\t%d\t%d\t%d\t%d\t%d\t%d\n", totalA, totalT, totalC, totalG, totalGap, totalN)
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"faInfo - Returns summary statistics to standard out for an\n" +
			"input.fa, including the counts for each base.\n" +
			"Usage:\n" +
			"faInfo inFile.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	var outFile *string = flag.String("outFile", "stdout", "redirect the output to a user specified file name, writes to stdout by default")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)

	faInfo(inFile, *outFile)
}
