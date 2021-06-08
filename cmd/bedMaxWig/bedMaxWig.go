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

func bedMaxWig(infile string, database string, chromsizeFile string, outfile string) {
	var records []*bed.Bed = bed.Read(infile)
	var wigData []wig.Wig = wig.Read(database)
	var sizes []chromInfo.ChromInfo = chromInfo.ReadToSlice(chromsizeFile)
	var outlist []*bed.Bed
	var currentBed *bed.Bed = records[0]
	var chromSlice []float64
	var i int

	for i = 0; i < len(sizes); i++ {
		chromSlice = WigChromToSlice(wigData, sizes[i].Size, sizes[i].Name)
		for k := 0; k < len(records); k++ {
			if records[k].Chrom == sizes[i].Name {
				currentBed = records[k]
				currentBed.Annotation = append(currentBed.Annotation, fmt.Sprintf("%f", bedRangeMax(chromSlice, records[k].ChromStart, records[k].ChromEnd)))
				outlist = append(outlist, currentBed)
			}
		}
	}
	bed.Write(outfile, outlist, 7)
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

	bedMaxWig(infile, database, chromsize, outfile)
}
