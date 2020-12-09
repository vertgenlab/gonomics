package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
)

func fqFormat(inFile string, outFile string) {
	records := fastq.Read(inFile)


	fastq.Write(outFile, records)
}

func usage() {
	fmt.Print(
		"fqFormat - reformat the sequences in a fastq file\n" +
			"Usage:\n" +
			" fqFormat input.fq output.fq\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	//var barcodeLength *int = flag.Int("barcodeLength", 16, "length of 10X cell barcode")
	//var umiLength *int = flag.Int("umiLength", 12, "length of transcript UMI")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	fqFormat(inFile, outFile)
}
