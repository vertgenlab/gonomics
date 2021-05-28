package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func vcfAfs(vcfFile string, outFile string, s popgen.McmcSettings) {
	g, err := popgen.VcfToAfs(vcfFile, s)
	exception.FatalOnErr(err)
	f := popgen.AfsToFrequency(*g)
	out := fileio.EasyCreate(outFile)
	defer out.Close()
	for i := 0; i < len(f); i++ {
		fmt.Fprintf(out, "%f\n", f[i])
	}
}

func usage() {
	fmt.Print(
		"vcfAfs - Returns allele frequency spectrum information in a text file for graphing.\n" +
			"vcfAfs - vcfFile.vcf outFile.txt\n")
}

func main() {
	var unPolarized *bool = flag.Bool("unPolarized", false, "vcfAfs creates polarized derived frequency spectra by default. When true, the cmd returns unpolarized site frequency spectra.")

	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	s := popgen.McmcSettings{
		UnPolarized: *unPolarized,
	}

	vcfFile := flag.Arg(0)
	outFile := flag.Arg(1)
	vcfAfs(vcfFile, outFile, s)
}
