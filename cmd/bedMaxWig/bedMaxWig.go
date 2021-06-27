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

func bedMaxWig(infile string, database string, chromsizeFile string, outfile string, norm *bool) {
	var records []bed.Bed = bed.Read(infile)
	var wigData []wig.Wig = wig.Read(database)
	var sizes []chromInfo.ChromInfo = chromInfo.ReadToSlice(chromsizeFile)
	var outlist []bed.Bed
	var currentBed bed.Bed = records[0]
	var chromSlice []float64
	var i int
	var wigTotal float64 = 0

	if *norm == true {
		for i := 0; i < len(wigData); i++ { //Goal here is to cycle through all the "chromosomes" of the wig: wig[0], wig[1], etc.
			var wigCounterByChrom float64 = 0
			for k := 0; k < len(wigData[i].Values); k++ { //Cycle through each value in the float64[]
				var chromValueMultiplyByStep float64 = 0
				chromValueMultiplyByStep = float64(wigData[i].Step) * wigData[i].Values[i] // multiply each value by the step for that chrom
				wigCounterByChrom = wigCounterByChrom + chromValueMultiplyByStep
			}
			wigTotal = wigTotal + wigCounterByChrom
		}
	}


	for i = 0; i < len(sizes); i++ {
		chromSlice = WigChromToSlice(wigData, sizes[i].Size, sizes[i].Name)
		for k := 0; k < len(records); k++ {
			if records[k].Chrom == sizes[i].Name {
				currentBed = records[k]
				if *norm == true {
					maxWig := bedRangeMax(chromSlice, records[k].ChromStart, records[k].ChromEnd)
					normMaxWig := maxWig/wigTotal
					currentBed.Annotation = append(currentBed.Annotation, fmt.Sprintf("%e", normMaxWig))
				} else {
					currentBed.Annotation = append(currentBed.Annotation, fmt.Sprintf("%f", bedRangeMax(chromSlice, records[k].ChromStart, records[k].ChromEnd)))
					outlist = append(outlist, currentBed)
				}
			}
		}
	}
	bed.Write(outfile, outlist)
}

func WigChromToSlice(w []wig.Wig, size int, chrom string) []float64 {
	output := make([]float64, size)
	for _, v := range w {
		if v.Chrom == chrom {
			output = v.Values
		}
	}
	return output
}

func sliceRangeAverage(w []float64, start int, end int) float64 {
	length := end - start
	var sum float64
	var i int

	for i = 0; i < length; i++ {
		sum = sum + w[i+start]
	}
	return (sum / float64(length))
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
		"bedMaxWig - Returns annotated bed with max wig score in bed entry range.\n" +
			"Usage:\n" +
			"bedMaxWig input.bed database.wig chrom.sizes output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var norm *bool = flag.Bool("normalize", false, "When true, will normalize the bedMaxWig output by dividing value by total wig hits.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	database := flag.Arg(1)
	chromsize := flag.Arg(2)
	outfile := flag.Arg(3)

	if *norm {
		bedMaxWig(infile, database, chromsize, outfile, norm)
		return
	}
	bedMaxWig(infile, database, chromsize, outfile, norm)

}
