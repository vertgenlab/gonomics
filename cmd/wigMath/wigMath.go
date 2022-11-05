package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func wigMath(s Settings) {
	var subtract []wig.Wig
	var i, j, chromIndex int
	records := wig.Read(s.InFile)
	multipleOptionCheck(s)
	if s.ElementWiseSubtract != "" {
		subtract = wig.Read(s.ElementWiseSubtract)
		for i = range records {
			chromIndex = getChromIndex(subtract, records[i].Chrom)
			for j = range records[i].Values {
				records[i].Values[j] -= subtract[chromIndex].Values[j]
			}
		}
	}
	if s.MovingAverageSmoothing > 1 {
		records = wig.SmoothSlice(records, s.MovingAverageSmoothing)
	}
	wig.Write(s.OutFile, records)
}

func multipleOptionCheck(s Settings) {
	var optionCount = 0
	if s.ElementWiseSubtract != "" {
		optionCount++
	}
	if s.MovingAverageSmoothing > 1 {
		optionCount++
	}
	if optionCount > 1 {
		log.Fatalf("wigMath can perform only one mathematical operation at a time. Rerun with a single option.")
	}
}

func getChromIndex(w []wig.Wig, chrom string) int {
	for i := range w {
		if w[i].Chrom == chrom {
			return i
		}
	}
	log.Fatalf("Error. Chromosome %v not found in wig.", chrom)
	return -1
}

func usage() {
	fmt.Print(
	"wigMath - Perform mathematical operations on wig format data.\n" +
		"Mathematical operations must be performed as single operations.\n" +
	"Usage:\n" +
	"wigMath in.wig out.wig" +
	"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile string
	OutFile string
	ElementWiseSubtract string
	MovingAverageSmoothing int
}

func main() {
	var expectedNumArgs = 2
	var elementWiseSubtract *string = flag.String("elementWiseSubtract", "", "Specify a second wig file to subtract (element-wise), from the first.")
	var movingAverageSmoothing *int = flag.Int("movingAverageSmoothing", 1, "Set to a number greater than 1 to perform moving average smoothing on input wig data.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	s := Settings {
		InFile: inFile,
		OutFile: outFile,
		ElementWiseSubtract: *elementWiseSubtract,
		MovingAverageSmoothing: *movingAverageSmoothing,
	}
	wigMath(s)
}