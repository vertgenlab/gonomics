package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"log"
	"os"
)

func usage() {
	fmt.Print(
		"sequelOverlap - A tool to find non/overlapping genomic regions\n\n" +
			"Gonomics Software\n" +
			"Author: Lowe Lab\thttp://www.vertgenlab.org\n\n" +
			"Source code: https://github.com/vertgenlab/gonomics\n\n" +
			"Version: 0.1.0\n\n" +
			"Usage:\n" +
			"  sequelOverlap [options] target query\n\n" +
			"Options:\n" +
			"  axt\t\taxt alignment as target select file\n" +
			"  bed\t\tbed regions as target select file\n" +
			//"  sam\t\tsam mapped coordinates as target select file\n"+
			"  vcf\t\tvcf genomic variance as target select file\n" +
			"  help\t\tview detailed help message specified option\n\n" +
			"Settings:\n" +
			"  -e, --extend\t\tadd to start and end of target coordinates, requires a chrom.sizes file\n"+
			"  -f, --format\t\tchoose format of final output file\n" +
			"  -n, --nonOverlap\treturn non-overlapping regions\n" +
			"  -t, --threads\t\tnumber of CPUs for Goroutine concurrency\n" +
			"  -o, --output\t\tfilename of final data\n\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	if len(os.Args) < 2 {
		flag.Usage()
	} else {
		switch os.Args[1] {
		case "axt":
		case "bed":
		case "vcf":
		default:
			flag.Usage()
		}
	}
	//var vcfArg *string = flag.String("v", "", )

	//log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	//flag.Parse()

	//file := flag.Arg(0)
}

func AxtBedCompare(axtfile string, bedfile string, output string, overlapSelect bool, target bool, query bool) {
	log.Printf("Reading axt...")
	blastz := axt.Read(axtfile)
	log.Printf("Reading bed...")
	peaks := bed.Read(bedfile)
	if overlapSelect {
		axt.FilterOverlap(blastz, peaks, output, target, query)
	} else {
		axt.NonOverlapFilter(blastz, peaks, output, target, query)
	}
}

/*
func axtToBed(axtfile string, bedfile string, output string) {
	log.Printf("Reading axt...")
	blastz := axt.Read(axtfile)
	log.Printf("Reading bed...")
	peaks := bed.Read(bedfile)
	axt.NonOverlapAxtBed(blastz, peaks, output)
}*/
