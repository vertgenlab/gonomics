// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func samToBed(samFilename string, bedFilename string, paired bool, fragLength int) {
	var aln sam.Sam
	var done bool = false

	//sam file to read
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)
	chroms := chromInfo.SliceToMap(header.Chroms)

	//bed file to write
	bedFile := fileio.EasyCreate(bedFilename)
	defer bedFile.Close()

	for aln, done = sam.ReadNext(samFile); done != true; aln, done = sam.ReadNext(samFile) {
		if aln.Cigar[0].Op != '*' {
			if fragLength != -1 {
				bed.WriteToFileHandle(bedFile, convert.SamToBedFrag(aln, fragLength, chroms), 4)
			} else {
				bed.WriteToFileHandle(bedFile, convert.SamToBed(aln), 4)
			}
		}
	}
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
	var paired *bool = flag.Bool("pairedEnd", false, "Specifies paired end reads")
	var fragLength *int = flag.Int("fragLength", -1, "Specifies the fragment length for ChIP-Seq")

	flag.Usage = usage
	flag.Parse()

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	log.Printf("Paired: %t\n", *paired)
	log.Printf("fragLength: %d\n", *fragLength)

	/*if paired && (*fragLength != -1) {
	        log.Fatalf("Error: cannot be both paired and have a fixed frag size.")
	}*/

	samToBed(inFile, outFile, *paired, *fragLength)
}
