package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math"
)

func vcfFilter(infile string, outfile string, chrom string, minPos int64, maxPos int64, ref string, alt string, minQual float64) {
	ch, _ := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	defer out.Close()

	for v := range ch {
		if vcf.Filter(v, chrom, minPos, maxPos, ref, alt, minQual) {
			vcf.WriteVcf(out.File, v)
		}
	}
}

func usage() {
	fmt.Print(
		"vcfFilter\n" +
			"Usage:\n" +
			"vcfFilter input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var minPos *int64 = flag.Int64("minPos", math.MinInt64, "Specifies the minimum position of the variant.")
	var maxPos *int64 = flag.Int64("maxPos", math.MaxInt64, "Specifies the maximum position of the variant.")
	var minQual *float64 = flag.Float64("minQual", 0.0, "Specifies the minimum quality score.")
	var ref *string = flag.String("ref", "", "Specifies the reference field.")
	var alt *string = flag.String("alt", "", "Specifies the alt field.")

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

	vcfFilter(infile, outfile, *chrom, *minPos, *maxPos, *ref, *alt, *minQual)
}
