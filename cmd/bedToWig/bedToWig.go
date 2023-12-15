// Command Group: "Data Conversion"

// Converts bed score to wig
package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/wig"
)

func bedToWig(method string, inFile string, refFile string, outFile string, missing float64, useRange bool, annotationField int) {
	ref := chromInfo.ReadToMap(refFile)
	var outWig []wig.Wig
	if method == "Reads" {
		rec := bed.Read(inFile)
		outWig = convert.BedReadsToWig(rec, ref)
	} else if method == "Name" || method == "Score" || method == "Annotation" {
		outWig = convert.BedValuesToWig(inFile, ref, missing, method, useRange, annotationField)
	} else {
		log.Fatalf("Unrecognized method. Expected 'Reads', 'Name', 'Score', or 'Annotation'. Found: %s.", method)
	}
	wig.SortByCoord(outWig)
	wig.Write(outFile, outWig)
}

func usage() {
	fmt.Print(
		"bedToWig - Converts bed score to wig\n" +
			"Usage:\n" +
			"bedToWig method input.bed reference.chrom.sizes output.wig\n" +
			"Method must be one of the following:\n" +
			"Score: Use the bed score column to set the wig value at the bed entry midpoint.\n" +
			"Reads: Use the bed region count to set the wig values across the entire range of the bed entry.\n" +
			"Name: Use the bed name column to set the wig value at the bed entry midpoint.\n" +
			"Annotation: Use an annotation column to set wig values at the bed entry midpoint. Default first annotation column, but can be controlled by the option 'annotationField'.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var missing *float64 = flag.Float64("missingData", 0, "For BedNameToWig, sets the value of the output wig in regions where there is no bed data.")
	var useRange *bool = flag.Bool("useRange", false, "For Name, Annotation, or Score method, set the wig value across the whole bed range instead of the midpoint.")
	var annotationField *int = flag.Int("annotationField", 0, "Specify which annotation column to use for wig values.")

	flag.Usage = usage

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	method := flag.Arg(0)
	inFile := flag.Arg(1)
	reference := flag.Arg(2)
	outFile := flag.Arg(3)

	bedToWig(method, inFile, reference, outFile, *missing, *useRange, *annotationField)
}
