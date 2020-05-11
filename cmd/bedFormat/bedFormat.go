package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func bedFormat(infile string, outfile string, fields int, ensemblToUCSC bool, UCSCToEnsembl bool) {
	ch := bed.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)

	if ensemblToUCSC && UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatable.")
	}

	for v := range ch {
		if ensemblToUCSC {
			v.Chrom = convert.EnsemblToUCSC(v.Chrom)
		}
		if UCSCToEnsembl {
			v.Chrom = convert.UCSCToEnsembl(v.Chrom)
		}
		bed.WriteBed(out.File, v, fields)
	}
}

func usage() {
	fmt.Print(
		"bedFormat: Options alter bed formatting.\n" +
			"Usage:\n" +
			"bedFormat input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")
	var fields *int = flag.Int("fields", 4, "Specify the number of bed fields.")

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

	bedFormat(infile, outfile, *fields, *ensemblToUCSC, *UCSCToEnsembl)
}
