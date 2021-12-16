// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func samToBed(samFilename string, bedFilename string, fragLength int) {
	var aln sam.Sam
	var done bool = false

	//sam file to read
	samFile := fileio.EasyOpen(samFilename)
	var err error

	header := sam.ReadHeader(samFile)
	chroms := chromInfo.SliceToMap(header.Chroms)

	//bed file to write
	bedFile := fileio.EasyCreate(bedFilename)

	for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
		if aln.Cigar[0].Op != '*' {
			if fragLength != -1 {
				bed.WriteToFileHandle(bedFile, convert.SamToBedFrag(aln, fragLength, chroms))
			} else {
				bed.WriteToFileHandle(bedFile, convert.SamToBed(aln))
			}
		}
	}
	err = samFile.Close()
	exception.PanicOnErr(err)
	err = bedFile.Close()
	exception.PanicOnErr(err)

	/* TODO: Write paired command
	if paired {
		outBed = convert.SamToBedPaired(records)
	} else*/
	/*
		if fragLength != int64(-1) {
			outBed = convert.SamToBedFrag(records, fragLength, ref)
		} else {
			outBed = convert.SamToBed(records)
		}*/
	/* When paired command is ready, add a flag in main
	var paired *bool = flag.Bool("pairedEnd", false, "Specifies paired end reads")
	log.Printf("Paired: %t\n", *paired)
	if paired && (*fragLength != -1) {
					log.Fatalf("Error: cannot be both paired and have a fixed frag size.")
	}
	samToBed(inFile, outFile, *paired, *fragLength)
	Change function samToBed to include paired flag
	Update samToBed_test.go
	*/
}

func usage() {
	fmt.Print(
		"samToBed - Converts sam to bed\n" +
			"Usage:\n" +
			" samToBed input.sam output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var fragLength *int = flag.Int("fragLength", -1, "Specifies the fragment length for ChIP-Seq, must be greater than or equal to read length")

	flag.Usage = usage
	flag.Parse()

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	log.Printf("fragLength: %d\n", *fragLength)

	samToBed(inFile, outFile, *fragLength)
}
