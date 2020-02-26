package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafFilter(inFile string, outFile string, threshold float64) {
	mafRecords := maf.Read(inFile)
	var outMaf []*maf.Maf

	for i, _ := range mafRecords {
		if mafRecords[i].Score >= threshold {
			outMaf = append(outMaf, mafRecords[i])
		}
	}

	maf.Write(outFile, outMaf)
}

func usage() {
	fmt.Print(
		"mafFilter - Filter a maf file to remove entries below a score threshold\n" +
			"Usage:\n" +
			" mafFilter mafFile oufMaf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	var threshold *float64 = flag.Float64("threshold", 0, "Specifies the threshold value")
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mafFile := flag.Arg(0)
	outMaf := flag.Arg(1)

	mafFilter(mafFile, outMaf, *threshold)
}
