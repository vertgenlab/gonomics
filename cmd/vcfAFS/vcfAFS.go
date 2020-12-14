package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func vcfAFS(vcfFile string, outFile string) {
	g := popgen.VcfToAFS(vcfFile)
	f := popgen.AFSToFrequency(g)
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	for i := 0; i < len(f); i++ {
		fmt.Fprintf(out, "%f\n", f[i])
	}
}

func usage() {
	fmt.Print(
		"vcfAFS - Returns allele frequency spectrum information in a text file for graphing.\n" +
			"vcfAFS - vcfFile.vcf outFile.txt\n")
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	vcfAFS(vcfFile, outFile)
}
