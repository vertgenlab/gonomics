package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func fastqFormat(s Settings) {
	if !s.PairedEnd {
		log.Fatalf("fastqFormat is still under development. Currently, the only formatting options available are for paired end reads. Select 'pairedEnd' from options.")
	} else if s.SingleCell {
		scPairChan := make(chan fastq.SingleCellPair)
		go fastq.ReadToChanSingleCellPair(s.R1InFile, s.R2InFile, s.BarcodeLength, s.UmiLength, scPairChan)
		outR1 := fileio.EasyCreate(s.R1OutFile)
		defer outR1.Close()
		outR2 := fileio.EasyCreate(s.R2OutFile)
		defer outR2.Close()

		for i := range scPairChan {
			fastq.WriteToFileHandle(outR1, i.Reads.Fwd)
			fastq.WriteToFileHandle(outR2, i.Reads.Rev)
		}
	} else {
		log.Fatalf("fastqFormat is still under development. Currently, the only formatting option available is to convert regular fastq reads into 10x formatted fastqs. Selected 'singleCell' from options.")
	}
}

func usage() {
	fmt.Print(
		"fastqFormat - Options alter fastq file formatting.\n" +
			"Usage:\n" +
			"fastqFormat input.fq output.fq\n" +
			"OR\n" +
			" fastqFormat -pairedEnd inputFwd.fq inputRev.fq outputFwd.fq outputRev.fq\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile        string
	OutFile       string
	R1InFile      string
	R2InFile      string
	R1OutFile     string
	R2OutFile     string
	PairedEnd     bool
	SingleCell    bool
	BarcodeLength int
	UmiLength     int
}

func main() {
	var expectedNumArgs int = 2
	var pairedEnd *bool = flag.Bool("pairedEnd", false, "Paired end reads, use two input and output fastq files.")
	var singleCell *bool = flag.Bool("singleCell", false, "Adds the single-cell barcode and UMI to the name of paired-end reads in 10X format.")
	var barcodeLength *int = flag.Int("barcodeLength", 16, "length of 10X cell barcode.")
	var umiLength *int = flag.Int("umiLength", 12, "length of transcript UMI.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *pairedEnd {
		expectedNumArgs = 4
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	s := Settings{
		InFile:        "",
		OutFile:       "",
		R1InFile:      "",
		R2InFile:      "",
		R1OutFile:     "",
		R2OutFile:     "",
		PairedEnd:     *pairedEnd,
		SingleCell:    *singleCell,
		BarcodeLength: *barcodeLength,
		UmiLength:     *umiLength,
	}

	if s.PairedEnd {
		s.R1InFile = flag.Arg(0)
		s.R2InFile = flag.Arg(1)
		s.R1OutFile = flag.Arg(2)
		s.R2OutFile = flag.Arg(3)
	} else {
		s.InFile = flag.Arg(0)
		s.OutFile = flag.Arg(1)
	}

	fastqFormat(s)
}
