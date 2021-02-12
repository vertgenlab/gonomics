package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
  "github.com/vertgenlab/gonomics/wig"
	"log"
)

func wigPeaks(in_wig string, out_bed string, threshold float64) {//threshold is float64 because WigValue Value aka v2.Value is float64. Need to cast int64 for bed.Bed for now

  //check Go documentation for gonomics > wig > functions written for wig
  records := wig.Read(in_wig) //type is []*Wig, aka slice of pointers to Wig
  inPeak := false //a bool to keep track of whether we are inside a peak
  var answer []*bed.Bed = make([]*bed.Bed,0) //to store list of peaks as slice of beds, need legnth argument start iwh 0
  var current bed.Bed //to store the current peak as a bed entry

  for _,v1 := range records { //in range for loop, i is index (0,1,2..) which we don't use in this instance, v is value (content of slice). Record is slice of wig struct, each iteration is slice/object of wig struct, which is often organized by chromosomes (or part of chromosome), and each wig struct will produce independent peaks. The v1 iteration has chr=chrom3,step=100, and a list of WigValues where there are positions+values, values like 11 22 100 etc. 
    for _,v2 := range v1.Values { //each v1 is a wig struct, whose Values is a []*WigValue which contains many Position-Value pairs. Each i2 is index which we don't use, and each v2 is a slice of WigValue, which contains .Position and .Value

      if v2.Value >= threshold { //either (from outside of a peak) start a new peak or inside of a peak
        if !inPeak { //this means start a new peak
          inPeak = true
          current = bed.Bed{Chrom: v1.Chrom, ChromStart: int64(v2.Position), ChromEnd: int64(v2.Position), Name: "", Score: int64(v2.Value)} //ChromEnd will correct once we reach the end of peak (in this program, the end of peak is still > threshold) OR end at >threshold in a peak, Score is the max Wig of the bed region, will update when inside the peak
        } else { //this means already inside peak
          current.ChromEnd = int64(v2.Position) //Update ChromEnd
          if int64(v2.Value) > current.Score { //Update Score if found new max wig score
            current.Score = int64(v2.Value)
          }
        }

      } else { //either (from inside of a peak) ending a peak or outside of a peak
        if inPeak { //if inside of a peak ending a peak
          inPeak = false
          answer = append(answer, &current) //pointer to current bed, add this pointer to answer which is also a pointer to bed slice (a collection of beds)
        }
        //if outside of a peak, nothing to do
      }

    }
  }

  bed.Write(out_bed, answer, 5)
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
  peakThreshold := flag.Float64("threshold",20,"if number of reads >= threshold, start calling peak")
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	in_wig := flag.Arg(0)
	out_bed := flag.Arg(1)
  threshold := *peakThreshold

	wigPeaks(in_wig, out_bed, threshold)
}
