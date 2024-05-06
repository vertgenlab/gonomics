package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"os"
	"strings"
)

// ToBedSettings defines the usage settings for the gtfTools toBed subcommand.
type ToBedSettings struct {
	InFile             string
	OutFile            string
	Tss                bool
	FirstTwoCodonBases bool
	ChromSizeFile      string
	Merge              bool
}

// toBedUsage defines the usage statement for the gtfTools toBed subcommand.
func toBedUsage(toBedFlags *flag.FlagSet) {
	fmt.Printf("gtfTools toBed - a tool to convert a gtf file into bed format.\n" +
		"Usage:\n" +
		"gtfTools toBed in.gtf out.bed\n" +
		"options:\n")
	toBedFlags.PrintDefaults()
}

// parseToBedArgs is the main function for the gtfTools toBed subcommand. It parses flags,
// checks for errors, and executes the ToBed function.
func parseToBedArgs() {
	var expectedNumArgs int = 2
	var err error
	toBedFlags := flag.NewFlagSet("toBed", flag.ExitOnError)
	var tss *bool = toBedFlags.Bool("tss", false, "Return a bed of tss positions annotated only with the geneName. Must provide chrom sizes file.")
	var firstTwoCodonBases *bool = toBedFlags.Bool("firstTwoCodonBases", false, "Return a bed of the first two positions of each coding exon.")
	var chromSizeFile *string = toBedFlags.String("chromSizeFile", "", "Specifies the name of a chrom.sizes file.")
	var merge *bool = toBedFlags.Bool("merge", false, "Merge overlapping entries after converting all records to beds. Available with tss only.")
	err = toBedFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	toBedFlags.Usage = func() { toBedUsage(toBedFlags) }
	if len(toBedFlags.Args()) != expectedNumArgs {
		toBedFlags.Usage()
		log.Fatalf("Error: expected %d arguments, but got %d.\n", expectedNumArgs, len(toBedFlags.Args()))
	}

	inFile := toBedFlags.Arg(0)
	outFile := toBedFlags.Arg(1)
	s := ToBedSettings{
		InFile:             inFile,
		OutFile:            outFile,
		Tss:                *tss,
		ChromSizeFile:      *chromSizeFile,
		Merge:              *merge,
		FirstTwoCodonBases: *firstTwoCodonBases,
	}

	toBed(s)
}

// toBed parses an input gtf format file and converts each transcript into a bed entry.
func toBed(s ToBedSettings) {
	var line string
	var doneReading bool
	var nameString string
	var currBed bed.Bed
	var err error

	if s.Tss && s.FirstTwoCodonBases {
		log.Fatalf("Error: user cannot specify both tss and firstTwoBasesOfCodons.")
	}

	if s.Tss && s.ChromSizeFile == "" {
		log.Fatalf("Error: user must specify a chromSizes file to convert to a Tss bed.\n")
	}

	if s.Tss {
		records := gtf.Read(s.InFile)
		sizes := chromInfo.ReadToMap(s.ChromSizeFile)
		beds := gtf.GenesToTssBed(records, sizes, s.Merge)
		bed.Write(s.OutFile, beds)
	} else if s.FirstTwoCodonBases {
		records := gtf.Read(s.InFile)
		beds := gtf.GenesToBedFirstTwoCodonBases(records)
		bed.Write(s.OutFile, beds)
	} else {
		file := fileio.EasyOpen(s.InFile)
		out := fileio.EasyCreate(s.OutFile)
		for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
			words := strings.Split(line, "\t")
			nameString = words[1] + ":" + words[2]
			for i := 5; i < len(words); i++ {
				nameString = nameString + ":" + words[i]
			}
			currBed = bed.Bed{Chrom: words[0], ChromStart: parse.StringToInt(words[3]) - 1, ChromEnd: parse.StringToInt(words[4]), Name: nameString, Score: 0, Strand: bed.Positive, FieldsInitialized: 6}
			if words[6] == "-" {
				currBed.Strand = bed.Negative
			}
			bed.WriteBed(out, currBed)
		}
		err = file.Close()
		exception.PanicOnErr(err)
		err = out.Close()
		exception.PanicOnErr(err)
	}
}
