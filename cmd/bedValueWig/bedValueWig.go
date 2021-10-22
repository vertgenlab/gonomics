// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func bedValueWig(s Settings) {
	if s.MinFlag && s.AverageFlag {
		log.Fatalf("Cannot select both min and average in the same operation.")
	}

	var records []bed.Bed = bed.Read(s.Infile)
	var wigData []wig.Wig = wig.Read(s.WigFile)
	var sizes []chromInfo.ChromInfo = chromInfo.ReadToSlice(s.SizesFile)
	var outList []bed.Bed
	var currentBed bed.Bed = records[0]
	var currValue float64
	var chromIndex int
	var i int
	var wigTotal float64 = 0

	if s.TrimLeft > 0 || s.TrimRight > 0 {
		bed.Trim(records, s.TrimLeft, s.TrimRight)
	}

	if s.NormFlag == true {
		for i := range wigData { //Goal here is to cycle through all the "chromosomes" of the wig: wig[0], wig[1], etc.
			var wigCounterByChrom float64 = 0
			for k := range wigData[i].Values { //Cycle through each value in the float64[]
				var chromValueMultiplyByStep float64 = 0
				chromValueMultiplyByStep = float64(wigData[i].Step) * wigData[i].Values[k] // multiply each value by the step for that chrom
				wigCounterByChrom = wigCounterByChrom + chromValueMultiplyByStep
			}
			wigTotal += wigCounterByChrom
		}
	}

	for i = 0; i < len(sizes); i++ {
		chromIndex = getWigIndex(wigData, sizes[i].Name)
		for k := 0; k < len(records); k++ {
			if records[k].Chrom == sizes[i].Name {
				currentBed = records[k]
				if currentBed.FieldsInitialized < 7 {
					currentBed.FieldsInitialized = 7
				}
				if s.MinFlag {
					currValue = bedRangeMin(wigData[chromIndex].Values, records[k].ChromStart, records[k].ChromEnd)
				} else if s.AverageFlag {
					currValue = bedRangeAverage(wigData[chromIndex].Values, records[k].ChromStart, records[k].ChromEnd)
				} else {
					currValue = bedRangeMax(wigData[chromIndex].Values, records[k].ChromStart, records[k].ChromEnd)
				}
				if s.NormFlag == true {
					currValue = currValue / wigTotal
				}
				currentBed.Annotation = append(currentBed.Annotation, fmt.Sprintf("%g", currValue)) // %g will
				// print %e for large exponents, %f otherwise.
				// Precision for %g; it is the smallest number of digits necessary to identify the value uniquely.
				outList = append(outList, currentBed)
			}
		}
	}
	bed.Write(s.OutFile, outList)
}

func getWigIndex(w []wig.Wig, chrom string) int {
	for i, v := range w {
		if v.Chrom == chrom {
			return i
		}
	}
	log.Fatalf("Chromosome with name: %v not found in wig.", chrom)
	return -1
}

func bedRangeAverage(w []float64, start int, end int) float64 {
	length := end - start
	var sum float64
	var i int

	for i = 0; i < length; i++ {
		sum = sum + w[i+start]
	}
	return sum / float64(length)
}

func bedRangeMin(w []float64, start int, end int) float64 {
	var min = w[start]
	for i := start; i < end; i++ {
		min = numbers.MinFloat64(min, w[i])
	}
	return min
}

func bedRangeMax(w []float64, start int, end int) float64 {
	var max float64
	for i := start; i < end; i++ {
		max = numbers.MaxFloat64(max, w[i])
	}
	return max
}

func usage() {
	fmt.Print(
		"bedValueWig - Returns bed file with entries annotated based on the values corresponding to the region in a wig file. Currently can only handle fixedStep, Start = 1, Step =1 wig files.\n" +
			"bedValueWig returns the maximum wig value overlapping a bed region by default.\n" +
			"Usage:\n" +
			"bedValueWig input.bed database.wig chrom.sizes output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	Infile      string
	WigFile     string
	SizesFile   string
	OutFile     string
	NormFlag    bool
	AverageFlag bool
	MinFlag     bool
	TrimLeft    int
	TrimRight   int
}

func main() {
	var expectedNumArgs int = 4
	var min *bool = flag.Bool("min", false, "Annotate bed entries with the minimum wig value instead of the maximum.")
	var average *bool = flag.Bool("average", false, "Annotate bed entries with the average wig value instead of the maximum.")
	var norm *bool = flag.Bool("normalize", false, "When true, will normalize the bedMaxWig output by dividing value by total wig hits.")
	var trimLeft *int = flag.Int("trimLeft", 0, "Exclude values on the left most edge of the bed region for consideration.")
	var trimRight *int = flag.Int("trimRight", 0, "Exclude values on the right most edge of the bed region for consideration.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	database := flag.Arg(1)
	chromSize := flag.Arg(2)
	outfile := flag.Arg(3)

	s := Settings{
		Infile:      infile,
		WigFile:     database,
		SizesFile:   chromSize,
		OutFile:     outfile,
		NormFlag:    *norm,
		AverageFlag: *average,
		MinFlag:     *min,
		TrimLeft:    *trimLeft,
		TrimRight:   *trimRight,
	}

	bedValueWig(s)
}
