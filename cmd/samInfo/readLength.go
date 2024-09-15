package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
)

type ReadLengthSettings struct {
	InFile  string
	OutFile string
}

func readLengthUsage(readLengthFlags *flag.FlagSet) {
	fmt.Printf("samInfo readLength - a tool to generate read length distributions from SAM/BAM files.\n" +
		"Usage:\n" +
		"samInfo readLength in.sam/bam out.txt" +
		"options:\n")
}

func parseReadLengthArgs() {
	var expectedNumArgs int = 2
	var err error
	readLengthFlags := flag.NewFlagSet("readLength", flag.ExitOnError)
	err = readLengthFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	readLengthFlags.Usage = func() { readLengthUsage(readLengthFlags) }
	if len(readLengthFlags.Args()) != expectedNumArgs {
		readLengthFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(readLengthFlags.Args()))
	}

	inFile := readLengthFlags.Arg(0)
	outFile := readLengthFlags.Arg(1)
	s := ReadLengthSettings{
		InFile:  inFile,
		OutFile: outFile,
	}
	readLength(s)
}

func readLength(s ReadLengthSettings) {
	var err error
	var currReadLength int
	outHist := fileio.EasyCreate(s.OutFile)
	data, _ := sam.GoReadToChan(s.InFile)
	histogram := make([]int, 200)

	for read := range data {
		currReadLength = len(read.Seq)
		if currReadLength > len(histogram) {
			histogramBuffer := make([]int, currReadLength+10)
			copy(histogramBuffer, histogram)
			histogram = histogramBuffer
		}
		histogram[currReadLength]++
	}
	_, err = fmt.Fprintf(outHist, "ReadLength\tCount\n")
	exception.PanicOnErr(err)
	for currReadLength = range histogram {
		_, err = fmt.Fprintf(outHist, "%v\t%v\n", currReadLength, histogram[currReadLength])
		exception.PanicOnErr(err)
	}

	err = outHist.Close()
	exception.PanicOnErr(err)
}
