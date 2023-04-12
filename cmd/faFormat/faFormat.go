// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

type Settings struct {
	InFile      string
	OutFile     string
	LineLength  int
	NamesFile   string
	TrimName    bool
	ToUpper     bool
	RevComp     bool
	NoGaps      bool
	NoGapBed    string
	Index       bool
	MaskInvalid bool
}

func faFormat(s Settings) {
	var records []fasta.Fasta
	if s.MaskInvalid {
		records = fasta.ReadForced(s.InFile)
	} else {
		records = fasta.Read(s.InFile)
	}
	var exist bool
	var names []string
	namesMap := make(map[string]int)

	if s.NoGapBed != "" {
		beds := bed.UngappedRegionsAllFromFa(records)
		bed.Write(s.NoGapBed, beds)
	}

	if s.NoGaps {
		fasta.RemoveGaps(records)
	}

	if s.NamesFile != "" {
		names = fileio.Read(s.NamesFile)
		for i := range names {
			namesMap[names[i]] = 1
		}
	}
	for i := range records {
		if s.NamesFile != "" {
			if _, exist = namesMap[records[i].Name]; !exist {
				continue
			}
		}
		if s.TrimName {
			records[i] = fasta.TrimName(records[i])
		}
		if s.ToUpper {
			fasta.ToUpper(records[i])

		}
		if s.RevComp {
			fasta.ReverseComplement(records[i])
			records[i].Name = records[i].Name + "_RevComp"
		}

	}

	file := fileio.EasyCreate(s.OutFile)
	fasta.WriteToFileHandle(file, records, s.LineLength)
	err := file.Close()
	exception.PanicOnErr(err)

	if s.Index {
		idx := fasta.CreateIndex(s.OutFile)
		idxFile := fileio.EasyCreate(s.OutFile + ".fai")
		_, err = fmt.Fprint(idxFile, idx)
		exception.PanicOnErr(err)
		err = idxFile.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print(
		"faFormat - reformat the sequences in a fasta file\n" +
			"Usage:\n" +
			" faFormat input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var lineLength *int = flag.Int("lineLength", 50, "wrap sequence lines after this many characters")
	var fastaNamesFile *string = flag.String("fastaNamesFile", "", "Text file, each line of the file are the fasta records to be manipulated (trimName, toUpper, revComp). Default is all fasta entries.")
	var trimName *bool = flag.Bool("trimName", false, "if a fasta name contains spaces, retains only the first space delimited field")
	var toUpper *bool = flag.Bool("toUpper", false, "Convert all DNA bases to upper case.")
	var revComp *bool = flag.Bool("revComp", false, "Return the reverse complement for each sequence.")
	var noGaps *bool = flag.Bool("noGaps", false, "Remove gaps from all input sequences.")
	var noGapBed *string = flag.String("noGapBed", "", "Find genomic coordinates containing regions outside gaps and write to a user-specified bed filename.")
	var createIndex *bool = flag.Bool("index", false, "Create index file (outputs to output.fa.fai).")
	var maskInvalid *bool = flag.Bool("maskInvalid", false, "N-mask extended IUPAC nucleotides (includes UWSMKRYBDHV).")

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
		InFile:      inFile,
		OutFile:     outFile,
		LineLength:  *lineLength,
		NamesFile:   *fastaNamesFile,
		TrimName:    *trimName,
		RevComp:     *revComp,
		ToUpper:     *toUpper,
		NoGaps:      *noGaps,
		NoGapBed:    *noGapBed,
		Index:       *createIndex,
		MaskInvalid: *maskInvalid,
	}

	faFormat(s)
}
