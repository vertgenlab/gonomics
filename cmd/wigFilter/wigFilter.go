package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

type Settings struct {
	InFile string
	OutFile string
	Chrom string
}

func wigFilter(s Settings) {
	var pass bool
	records := wig.Read(s.InFile)
	var answer []wig.Wig = make([]wig.Wig, 0)

	for i := range records {
		pass = true
		if s.Chrom != "" && records[i].Chrom != s.Chrom {
			pass = false
		}
		if pass {
			answer = append(answer, records[i])
		}
	}
	wig.Write(s.OutFile, answer)
}

func usage() {
	fmt.Print(
		"wigFilter - Returns a filtered wig based on option parameters.\n" +
			"Usage:\n" +
			"wigFilter input.wig output.wig\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var chrom *string = flag.String("chrom", "", "Retains wig entries with chromosome name matching this value.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	var s Settings = Settings {
		InFile: inFile,
		OutFile: outFile,
		Chrom: *chrom,
	}

	wigFilter(s)
}
