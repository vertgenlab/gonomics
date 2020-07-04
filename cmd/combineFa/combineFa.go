package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func usage() {
	fmt.Print(
		"combineFa - concatenate two or more fasta files into a single file\n" +
			"Usage:\n" +
			" combineFa -out file.fa input.fa input.fa...\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	//must provide at least 2 fasta files
	var expectedNumArgs int = 1
	//Using option out to be able to use *.fa in command line
	var outMerge *string = flag.String("out", "/dev/stdout", "output filename for merging `.fa`sta records")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if len(flag.Args()) >= expectedNumArgs {
		combineFastaFiles(flag.Args(), *outMerge)
	} else {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

}

func combineFastaFiles(files []string, outputFile string) {
	ans := fasta.NewPipe()
	ans.Wg.Add(1)
	go fasta.ReadMultiFilesToChan(ans, files)
	go fasta.WritingChannel(outputFile, ans.StdOut, ans.Wg)
	ans.Wg.Wait()
}
