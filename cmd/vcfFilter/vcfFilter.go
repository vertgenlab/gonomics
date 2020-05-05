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
	ch := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	defer out.Close()
	var pass bool = false

	for v := range ch {
		pass = true
		if v.Pos < minPos {
			pass = false
		}
		if v.Pos > maxPos {
			pass = false
		}
		if chrom != "" && v.Chr != chrom {
			pass = false
		}
		if ref != "" && v.Ref != ref {
			pass = false
		}
		if alt != "" && v.Alt != alt {
			pass = false
		}
		if v.Qual < minQual {
			pass = false
		}
		if pass {
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
