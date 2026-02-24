package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"path/filepath"
)

type mapqSettings struct {
	InFile     string
	OutFile    string
	BedRegions string
}

func mapqUsage(mapqFlags *flag.FlagSet) {
	fmt.Printf("samInfo mapq - a tool to generate mapping quality statistics from SAM/BAM files.\n" +
		"Unampped reads will be ignored.\n" +
		"Usage:\n" +
		"samInfo mapq [options] in.sam/bam histogram.txt\n" +
		"options:\n")
	mapqFlags.PrintDefaults()
}

func parseMapqArgs() {
	var expectedNumArgs int = 2
	var err error
	mapqFlags := flag.NewFlagSet("mapq", flag.ExitOnError)

	var bedfile *string = mapqFlags.String("bedfile", "", "Provide a file with bed regions. Only reads that overlap a bed regions will be analyzed. Only compatible with"+
		"bam files. Any alignment with any degree of overlap will be considered.")

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

	if *bedfile != "" && filepath.Ext(inFile) != ".bam" {
		//since the bed region option uses the bamseeker function, we check to see if the user input is a bam (not a sam)
		log.Fatalf("Error: -bedfile must be used with a bam alignment file input. Alignment file provided: %s\n", inFile)
	}

	s := mapqSettings{
		InFile:     inFile,
		OutFile:    outFile,
		BedRegions: *bedfile,
	}
	mapq(s)
}

func mapq(s mapqSettings) {
	var otherMapq []uint8
	hist := make([]int, 61)
	out := fileio.EasyCreate(s.OutFile)

	if s.BedRegions != "" {
		hist, otherMapq = bedRegionsOnly(s, hist, otherMapq)
	} else {
		hist, otherMapq = wholeGenome(s, hist, otherMapq)
	}

	writeHist(out, hist, otherMapq)
	exception.PanicOnErr(out.Close())
}

func wholeGenome(s mapqSettings, hist []int, otherMapq []uint8) ([]int, []uint8) {
	samfile, _ := sam.GoReadToChan(s.InFile)

	//loop through the alignment file
	for i := range samfile {
		hist, otherMapq = addToHist(i, hist, otherMapq)
	}
	return hist, otherMapq
}

func bedRegionsOnly(s mapqSettings, hist []int, otherMapq []uint8) ([]int, []uint8) {
	var j int
	var toAnalyze []sam.Sam
	//get the bed regions for analysis
	beds := bed.Read(s.BedRegions)

	//open up the alinment file with a bam reader (more efficient search for overlapped regions)
	br, _ := sam.OpenBam(s.InFile)
	bai := sam.ReadBai(s.InFile + ".bai")

	//loop through bed file
	for i := range beds {
		//for each region, get the overlapping alignments
		toAnalyze = sam.SeekBamRegion(br, bai, beds[i].Chrom, uint32(beds[i].ChromStart), uint32(beds[i].ChromEnd))
		for j = range toAnalyze {
			//add the overlapping alignments to histogram
			hist, otherMapq = addToHist(toAnalyze[j], hist, otherMapq)
		}
	}
	return hist, otherMapq
}

func addToHist(aln sam.Sam, hist []int, otherMapq []uint8) ([]int, []uint8) {
	var found bool
	var j int

	if sam.IsUnmapped(aln) {
		//ignore unmapped reads
		return hist, otherMapq
	}
	//check for unusual mapq score. For example, sometimes cellranger uses mapq 255
	if aln.MapQ > 60 || aln.MapQ < 0 {
		found = false
		for j = range otherMapq { //check to see if the unusual mapq has been seen before
			if aln.MapQ == otherMapq[j] {
				//if it has, add it to the histogram and break
				hist[61+j]++
				found = true
				break
			}
		}
		if !found {
			//if it hasn't, add to both the histogram and list of unusual MapQ's
			otherMapq = append(otherMapq, aln.MapQ)
			hist = append(hist, 1)
		}
	} else {
		//normal case for mapq 0-60, iterate the slice with the index corresponding to mapq
		hist[aln.MapQ]++
	}
	return hist, otherMapq
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
