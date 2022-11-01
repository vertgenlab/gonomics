// Command Group: Deep Learning

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faToPredictSet(s Settings) {
	var err error
	var i, j int
	records := fasta.Read(s.InFile)
	out := fileio.EasyCreate(s.OutFile)
	var currFa fasta.Fasta
	var lineToWrite string

	for i = range records {
		for j = 0; j < len(records[i].Seq)-s.WindowSize; j += s.Stride {
			currFa = fasta.Extract(records[i], j, j+s.WindowSize, fmt.Sprintf("%s:%v-%v", records[i].Name, j, j+s.WindowSize))
			dna.AllToUpper(currFa.Seq)
			lineToWrite = fmt.Sprintf("%s\t%s\n", currFa.Name, dna.BasesToString(currFa.Seq))
			_, err = fmt.Fprintf(out, lineToWrite)
			exception.PanicOnErr(err)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"faToPredictSet - Make deep learning prediction TSV files from input fasta format data.\n" +
			"Usage:\n" +
			" faToPredictSet input.fa output.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile     string
	OutFile    string
	WindowSize int
	Stride     int
}

func main() {
	var expectedNumArgs int = 2
	var windowSize *int = flag.Int("windowSize", 400, "Set the size of the sequence window.")
	var stride *int = flag.Int("stride", 1, "Sets the size of the stride between prediction windows.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	s := Settings{
		InFile:     inFile,
		OutFile:    outFile,
		WindowSize: *windowSize,
		Stride:     *stride,
	}

	faToPredictSet(s)
}
