// Command Group: "FASTA and Multi-FASTA Tools"

// Reformat the sequences in a fasta file
package main

import (
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

type Settings struct {
	InFile          string
	OutFile         string
	LineLength      int
	NamesFile       string
	TrimName        bool
	ToUpper         bool
	RevComp         bool
	NoGaps          bool
	NoGapBed        string
	Index           bool
	MaskInvalid     bool
	MultiFaNoGapBed string
	QuerySeqName    string
	ChromName       string
	Rename          string
}

func faFormat(s Settings) {
	var records []fasta.Fasta
	var words []string
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

	if s.MultiFaNoGapBed != "" {
		if s.QuerySeqName == "" {
			log.Fatalf("Error: to use multiFaNoGapBed, must specify querySeqName.\n")
		}
		if s.ChromName == "" {
			log.Fatalf("Error: to use multiFaNoGapBed, must specify chromName.\n")
		}
		beds := bed.MultiFaUngappedRegions(records, s.ChromName, s.QuerySeqName)
		bed.Write(s.MultiFaNoGapBed, beds)
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

	if s.Rename != "" {
		words = strings.Split(s.Rename, ",")
		if len(words) != 2 {
			log.Fatalf("Error: expected two fields, comma delimited, in -rename. Found: %v.\n", s.Rename)
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
		if s.Rename != "" {
			if records[i].Name == words[0] {
				records[i].Name = words[1]
			}
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
	var multiFaNoGapBed *string = flag.String("multiFaNoGapBed", "", "Find genomic coordinates containing regions \n"+
		"outside gaps and write to a user-specified bed filename. Ungapped regions are reported for an aligned query sequence whose name must be specified by\n"+
		"the option querySeqName. User must also specify a chromName.")
	var querySeqName *string = flag.String("querySeqName", "", "Specify the name of the sequence in the multiFa from which to generate a noGap file with multiFaNoGapBed.")
	var chromName *string = flag.String("chromName", "", "Specify the name of the chromosome in the multiFa for multiFaNoGapBed.")
	var createIndex *bool = flag.Bool("index", false, "Create index file (outputs to output.fa.fai).")
	var maskInvalid *bool = flag.Bool("maskInvalid", false, "N-mask extended IUPAC nucleotides (includes UWSMKRYBDHV).")
	var rename *string = flag.String("rename", "", "Rename a name field using comma delimited argument (ex. 'old,new'). Only one name field can be changed at a time.")

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
		InFile:          inFile,
		OutFile:         outFile,
		LineLength:      *lineLength,
		NamesFile:       *fastaNamesFile,
		TrimName:        *trimName,
		RevComp:         *revComp,
		ToUpper:         *toUpper,
		NoGaps:          *noGaps,
		NoGapBed:        *noGapBed,
		Index:           *createIndex,
		MaskInvalid:     *maskInvalid,
		MultiFaNoGapBed: *multiFaNoGapBed,
		QuerySeqName:    *querySeqName,
		ChromName:       *chromName,
		Rename:          *rename,
	}

	faFormat(s)
}
