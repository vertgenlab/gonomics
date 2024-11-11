package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func filterScBam(inbam, barcodes, outbam string) {
	var found bool
	var err error
	var bx interface{}

	aln, head := sam.GoReadToChan(inbam)

	out := fileio.EasyCreate(outbam)
	bw := sam.NewBamWriter(out, head)

	mp := createBxMap(barcodes)

	for i := range aln {

		bx, _, err = sam.QueryTag(i, "CB")
		exception.PanicOnErr(err)

		_, found = mp[bx.(string)]

		if found {
			sam.WriteToBamFileHandle(bw, i, 0)
		}
	}
	err = bw.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func createBxMap(infile string) map[string]int {
	mp := make(map[string]int)
	in := fileio.Read(infile)

	for i := range in {
		mp[in[i]] = 0
	}
	return mp
}

func usage() {
	fmt.Print("filterScBam -- filter a cellranger output alignment file based on a list of cell barcodes\n" +
		"Usage:\n" +
		"filterScBam in.bam barcodes.txt out.bam\n")
}

func main() {
	var expectedNumArgs int = 3
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.PrintDefaults()
		usage()
		log.Fatalf("Expecting 3 arguments, got %d", len(flag.Args()))
	}

	inbam := flag.Arg(0)
	bx := flag.Arg(1)
	outBam := flag.Arg(2)

	filterScBam(inbam, bx, outBam)

}
