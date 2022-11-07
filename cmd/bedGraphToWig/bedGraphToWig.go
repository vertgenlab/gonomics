package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func bedGraphToWig(inFile string, chromFile string, outFile string, missing float64) {
	ref := chromInfo.ReadToMap(chromFile)
	outWig := convert.BedGraphToWig(inFile, ref, missing)
	wig.Write(outFile, outWig)
}

func usage() {
	fmt.Print(
		"bedGraphToWig - Converts bedGraph to wig. Wig scores will be equal to the bedGraph dataValue field across the range of the bedGraph entry.\n" +
			"Usage:\n" +
			"bedGraphToWig input.bedGraph reference.chrom.sizes output.wig\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var missing *float64 = flag.Float64("missingData", 0, "Sets the value of the output wig in regions where there is no bedGraph data.")

	flag.Usage = usage

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	reference := flag.Arg(1)
	outFile := flag.Arg(2)

	bedGraphToWig(inFile, reference, outFile, *missing)
}
