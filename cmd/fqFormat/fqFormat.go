package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func fqFormat(inFileFwd string, inFileRev string, barcodeLength int, umiLength int, outFileFwd string, outFileRev string) {
	TenXPairChan := make(chan *fastq.TenXPair)
	go fastq.ReadToChanTenXPair(inFileFwd, inFileRev, barcodeLength, umiLength, TenXPairChan)

	fileFwd := fileio.EasyCreate(outFileFwd)
	defer fileFwd.Close()
	fileRev := fileio.EasyCreate(outFileRev)
	defer fileRev.Close()

	UMImap := make(map[string]int)

	var currentUmi string

	for i := range TenXPairChan {
		currentUmi = dna.BasesToString(i.Umi)
		count, found := UMImap[currentUmi]

		if !found {
			UMImap[currentUmi] = 1

			fastq.WriteToFileHandle(fileFwd, i.Fwd)
			fastq.WriteToFileHandle(fileRev, i.Rev)
		} else {
			UMImap[currentUmi] = count + 1
		}
	}

}

func usage() {
	fmt.Print(
		"fqFormat - reformat the sequences in a fastq file\n" +
			"Usage:\n" +
			" fqFormat inputFwd.fq inputRev.fq outputFwd.fq outputRev.fq\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var barcodeLength *int = flag.Int("barcodeLength", 16, "length of 10X cell barcode")
	var umiLength *int = flag.Int("umiLength", 12, "length of transcript UMI")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFileFwd := flag.Arg(0)
	inFileRev := flag.Arg(1)
	outFileFwd := flag.Arg(2)
	outFileRev := flag.Arg(3)

	fqFormat(inFileFwd, inFileRev, *barcodeLength, *umiLength, outFileFwd, outFileRev)
}
