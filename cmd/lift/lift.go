package main

import (
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/fileio"
)

func lift(chainFile string, inFile string, outFile string) {
	//first task, make tree from chainFile
	chainChan, chainHeader := chain.GoReadToChan(chainFile)
	var chainIntervals []interval.Interval
	for val := range selectChan {
		chainIntervals = append(chainIntervals, val)
	}
	tree := interval.BuildTree(chainIntervals)

	//second task, read in intervals, find chain, and convert to new interval
	inChan := interval.GoReadToLiftChan(inFile)

	outChan := HELPERFUNCTION()
	out :=  fileio.EasyCreate(outFile)
	WRITETOFILE(out, outChan)
}

func usage() {
	fmt.Print(
		"lift - Moves an interval interface compatable file format between assemblies.\n" +
			"Usage:\n" +
			"lift lift.chain inFile outFile\n" +
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

	chainFile := flag.Arg(0)
	inFile := flag.Arg(1)
	outFile := flag.Arg(2)
	lift(chain, inFile, outFile)
}