package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

func bedShear(inFile string, outFile string, fragmentSize int) {
	var err error
	var currInBed bed.Bed
	var currOutBed bed.Bed = bed.Bed{}
	var currStart int

	if fragmentSize < 1 {
		log.Fatalf("Error: fragmentSize must be a positive integer. Found: %v.\n", fragmentSize)
	}

	ch := bed.GoReadToChan(inFile)
	out := fileio.EasyCreate(outFile)

	for currInBed = range ch {
		currOutBed.Chrom = currInBed.Chrom
		currOutBed.ChromStart = currInBed.ChromStart
		currOutBed.ChromEnd = currInBed.ChromStart + fragmentSize
		currOutBed.FieldsInitialized = currInBed.FieldsInitialized
		if currInBed.FieldsInitialized > 3 {
			currOutBed.Name = currInBed.Name
		}
		if currInBed.FieldsInitialized > 4 {
			currOutBed.Score = currInBed.Score
		}
		if currInBed.FieldsInitialized > 5 {
			currOutBed.Strand = currInBed.Strand
		}
		if currInBed.FieldsInitialized > 6 {
			currOutBed.Annotation = currInBed.Annotation
		}

		for currStart = currInBed.ChromStart; currStart < currInBed.ChromEnd; currStart += fragmentSize {
			currOutBed.ChromStart = currStart
			currOutBed.ChromEnd = numbers.Min(currStart+fragmentSize, currInBed.ChromEnd)
			bed.WriteBed(out, currOutBed)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"bedShear - Split bed entries into smaller fragment bed entries.\n" +
			"Usage:\n" +
			"bedShear input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var fragmentSize *int = flag.Int("fragmentSize", 1, "Set the maximum size of output bed fragments.")

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

	bedShear(inFile, outFile, *fragmentSize)
}
