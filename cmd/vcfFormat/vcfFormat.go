package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfFormat(infile string, outfile string, ensemblToUCSC bool, UCSCToEnsembl bool) {
	ch := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	defer out.Close()

	if ensemblToUCSC && UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatable.")
	}

	for v := range ch {
		if ensemblToUCSC {
			v.Chr = convert.EnsemblToUCSC(v.Chr)
		}
		if UCSCToEnsembl {
			v.Chr = convert.UCSCToEnsembl(v.Chr)
		}
		vcf.WriteVcf(out.File, v)
	}
}

func usage() {
	fmt.Print(
		"vcfFormat: Options alter VCF formatting.\n" +
			"Usage:\n" +
			"vcfFormat input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	vcfFormat(infile, outfile, *ensemblToUCSC, *UCSCToEnsembl)
}
