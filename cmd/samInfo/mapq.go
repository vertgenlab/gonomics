package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
)

type mapqSettings struct {
	InFile  string
	OutFile string
}

func mapqUsage(readLengthFlags *flag.FlagSet) {
	fmt.Printf("samInfo mapq - a tool to generate mapping quality statistics from SAM/BAM files.\n" +
		"Unampped reads will be ignored.\n" +
		"Usage:\n" +
		"samInfo mapq in.sam/bam histogram.txt\n" +
		"options:\n")
}

func parseMapqArgs() {
	var expectedNumArgs int = 2
	var err error
	mapqFlags := flag.NewFlagSet("mapq", flag.ExitOnError)
	err = mapqFlags.Parse(os.Args[2:])

	exception.PanicOnErr(err)
	mapqFlags.Usage = func() { mapqUsage(mapqFlags) }
	if len(mapqFlags.Args()) != expectedNumArgs {
		mapqFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(mapqFlags.Args()))
	}

	inFile := mapqFlags.Arg(0)
	outFile := mapqFlags.Arg(1)
	s := mapqSettings{
		InFile:  inFile,
		OutFile: outFile,
	}
	mapq(s)
}

func mapq(s mapqSettings) {
	//mp := make(map[uint8]int)
	var otherMapq []uint8
	var j int
	var found bool
	hist := make([]int, 61)
	out := fileio.EasyCreate(s.OutFile)
	aln, _ := sam.GoReadToChan(s.InFile)

	//loop through the alignment file
	for i := range aln {
		if sam.IsUnmapped(i) {
			//ignore unmapped reads
			continue
		}
		//check for unusual mapq score. For example, sometimes cellranger uses mapq 255
		if i.MapQ > 60 || i.MapQ < 0 {
			found = false
			for j = range otherMapq { //check to see if the unusual mapq has been seen before
				if i.MapQ == otherMapq[j] {
					//if it has, add it to the histogram and break
					hist[61+j]++
					found = true
					break
				}
			}
			if !found {
				//if it hasn't, add to both the histogram and list of unusual MapQ's
				otherMapq = append(otherMapq, i.MapQ)
				hist = append(hist, 1)
			}
		} else {
			//normal case for mapq 0-60, iterate the slice with the index corresponding to mapq
			hist[i.MapQ]++
		}
	}
	writeHist(out, hist, otherMapq)
	exception.PanicOnErr(out.Close())
}

func writeHist(out *fileio.EasyWriter, hist []int, otherMapq []uint8) {
	fileio.WriteToFileHandle(out, "mapQ\tcount")
	for i := range hist {
		if i < 61 {
			fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%d", i, hist[i]))
		} else {
			fileio.WriteToFileHandle(out, fmt.Sprintf("%d\t%d", otherMapq[i-61], hist[i]))
		}
	}
}
