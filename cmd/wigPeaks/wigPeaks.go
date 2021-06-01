package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/wig"
	"log"
)

func wigPeaks(in_wig string, out_bed string, threshold float64) { //threshold is float64 because WigValue Value aka v2.Value is float64. Need to cast int64 for bed.Bed for now

	//check Go documentation for gonomics > wig > functions written for wig
	records := wig.Read(in_wig) //type is []Wig, aka slice of Wig structs
	inPeak := false             //a bool to keep track of whether we are inside a peak
	//var answer []*bed.Bed = make([]*bed.Bed,0) //to store list of peaks as slice of beds, use pointer to pass by reference and edit slice of beds, "make" command needs legnth argument so start with 0
	var current *bed.Bed              //to store the current peak as a bed entry
	out := fileio.EasyCreate(out_bed) //instead of answer, use "chan" method to write as we go. out has type EasyCreate defined as a gonomics struct, one of the fields is File which has type *os.File, needed for bed.WriteBed
	defer out.Close()

	for _, v1 := range records { //in range for loop, i is index (0,1,2..) which we don't use in this instance, v is value (content of slice). Record is slice of wig struct, each iteration is slice/object of wig struct, which looks like a block and is often organized by chromosomes (or part of chromosome), and each wig struct will produce independent peaks. The v1 iteration has chr=chrom3,step=100, and a list of WigValues where there are positions+values, values like 11 22 100 etc.
		inPeak = false                 //when entering a new wig struct, set inPeak to false
		var wigPosition = v1.Start
		for _, v2 := range v1.Values { //each v1 is a wig struct, whose Values is a []float64 which contains the value at each position. Each i2 is index which we don't use, and each v2 is a float64.
			if v2 >= threshold { //either (from outside of a peak) start a new peak or inside of a peak
				if !inPeak { //this means start a new peak if this is currently set to false.
					inPeak = true                                                                                                           //must remember to set inPeak to reflect Peak status
					current = &bed.Bed{Chrom: v1.Chrom, ChromStart: wigPosition, ChromEnd: wigPosition + 1, Name: "", Score: int(v2)} //ChromEnd is +1 because bed has [open,closed) interval (note: in this program, the position that ends the peak is still >threshold, not the first position that dips <threshold), Score is the max Wig of the bed region, will update when inside the peak
				} else { //this means already inside peak
					current.ChromEnd = wigPosition + 1     //Update ChromEnd
					if v2 > float64(current.Score) { //Update Score if found new max wig score
						current.Score = int(v2)
					}
				}
			} else { //either (from inside of a peak) ending a peak or outside of a peak
				if inPeak { //if inside of a peak ending a peak
					inPeak = false //must remember to set inPeak to reflect Peak status
					//answer = append(answer, current) //pointer to current bed, add this pointer to answer which is also a pointer to bed slice (a collection of beds)
					bed.WriteBed(out.File, current, 5) //instead of answer, use "chan" method to write as we go
				}
				//if outside of a peak, nothing to do
			}
			wigPosition = wigPosition + v1.Step
		}
		if inPeak { //after wig struct ends, if inPeak is still true (i.e. ended on a value>threshold), should end peak and append current
			inPeak = false //redundant since this will happen when enter a new wig struct
			//answer = append(answer, current)
			bed.WriteBed(out.File, current, 5) //instead of answer, use "chan" method to write as we go
		}
	}

	//bed.Write(out_bed, answer, 5) //use answer (the pointer to bed slice) to write bed file
}

func usage() {
	fmt.Print(
		"wigPeaks - takes wig file and finds peaks\n" +
			"Usage:\n" +
			" wigPeaks in.wig out.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	peakThreshold := flag.Float64("threshold", 20, "if number of reads >= threshold, start calling peak")
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	in_wig := flag.Arg(0)
	out_bed := flag.Arg(1)
	threshold := peakThreshold

	wigPeaks(in_wig, out_bed, *threshold)
}
