package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafToBed(mafFile string, outBed string, reference string) {
	mafRecords := maf.Read(mafFile)
	var bedList []*bed.Bed

	for i, _ := range mafRecords {
		for k, _ := range mafRecords[i].Species {
			assembly, chrom := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src)
			if assembly == reference {
				if mafRecords[i].Species[k].SLine != nil {
					current := bed.Bed{Chrom: chrom, ChromStart: mafRecords[i].Species[k].SLine.Start, ChromEnd: mafRecords[i].Species[k].SLine.Start + mafRecords[i].Species[k].SLine.Size, Name: "blank", Score: int(mafRecords[i].Score)}
					bedList = append(bedList, &current)
				}
			}
		}
	}
	bed.Write(outBed, bedList, 5)
}

func usage() {
	fmt.Print(
		"mafToBed - convert a maf alignment into a bed, where the bed score is the alignment score\n" +
			"Usage:\n" +
			" mafToBed mafFile outBed referenceSpeciesName\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mafFile := flag.Arg(0)
	outBed := flag.Arg(1)
	reference := flag.Arg(2)

	mafToBed(mafFile, outBed, reference)
}
