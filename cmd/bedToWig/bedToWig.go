// Command Group: "Data Conversion"

// Converts bed score to wig
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
)

type Settings struct {
	Method          string
	InFile          string
	RefFile         string
	OutFile         string
	DefaultValue    float64
	UseRange        bool
	AnnotationField int
}

func bedToWig(s Settings) {
	ref := chromInfo.ReadToMap(s.RefFile)
	var outWig = wig.MakeSkeleton(ref, s.DefaultValue)
	if s.Method == "Reads" {
		rec := bed.Read(s.InFile)
		outWig = convert.BedReadsToWig(rec, ref)
	} else if s.Method == "Name" || s.Method == "Score" || s.Method == "Annotation" {
		outWig = convert.BedValuesToWig(s.InFile, ref, s.DefaultValue, s.Method, s.UseRange, s.AnnotationField)
	} else {
		log.Fatalf("Unrecognized method. Expected 'Reads', 'Name', 'Score', or 'Annotation'. Found: %s.", s.Method)
	}
	wig.Write(s.OutFile, outWig)
}

func usage() {
	fmt.Print(
		"bedToWig - Converts bed score to wig\n" +
			"Usage:\n" +
			"bedToWig method input.bed reference.chrom.sizes output.wig\n" +
			"Method must be one of the following:\n" +
			"\tScore: Use the bed score column to set the wig value at the bed entry midpoint.\n" +
			"\tReads: Use the bed region count to set the wig values across the entire range of the bed entry.\n" +
			"\tName: Use the bed name column to set the wig value at the bed entry midpoint.\n" +
			"\tAnnotation: Use an annotation column to set wig values at the bed entry midpoint. Default first annotation column, but can be controlled by the option 'annotationField'.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var defaultValue *float64 = flag.Float64("defaultValue", math.MaxFloat64, "(Developer) Set the default value of the wig struct. Set to a value that is unlikely to be found in the real wig data.")
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

	s := Settings{
		Method:          method,
		InFile:          inFile,
		RefFile:         reference,
		OutFile:         outFile,
		DefaultValue:    *defaultValue,
		UseRange:        *useRange,
		AnnotationField: *annotationField,
	}

	bedToWig(s)
}
