package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/bed"
	"log"
)

func humanInsertion(mafFile string, outBed_ins string, outBed_del string, species_ins string, species_del string) {
	//initialize variables
	mafRecords := maf.Read(mafFile) //Read entire mafFile. mafRecords has type Maf
	var bedList_ins []*bed.Bed //initialize 2 bed files
	var bedList_del []*bed.Bed //1 bed file for ins, 1 bed file for del

	//go through each line
	for i, _ := range mafRecords {//each i is a block
		for k, _ := range mafRecords[i].Species {//each k is a line

			//convert maf to bed, start with getting assembly because it is needed to verify species_ins and species_del
			assembly_del, chrom_del := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src) //start with species_del, get assembly (e.g. panTro6), chrom (e.g. chrI)
			assembly_ins, chrom_ins := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k-1].Src) //then find corresponding species_ins line, same block, 1 line above since pairwise maf

			//get eC species_del lines
			if mafRecords[i].Species[k].ELine != nil {
				if mafRecords[i].Species[k].ELine.Status=='C' && assembly_del == species_del { //I decided to check for both Status and Src here because they are on the same data level

					//get corresponding s species_ins lines
					if assembly_ins != species_ins { //verify line k-1 is indeed species_ins
						log.Fatalf("species_ins was incorrect. Please check you have a pairwise maf file, and entered species_ins and species_del correctly") //otherwise fatal
					}
					if mafRecords[i].Species[k-1].SLine != nil { //if corresponding species_ins line an s line

						//convert maf to bed, continued
						current_del := bed.Bed{Chrom: chrom_del, ChromStart: mafRecords[i].Species[k].ELine.Start, ChromEnd: mafRecords[i].Species[k].ELine.Start + mafRecords[i].Species[k].ELine.Size, Name: "del", Score: int64(mafRecords[i].Score)} //get chrom,start,end,name,score
						current_ins := bed.Bed{Chrom: chrom_ins, ChromStart: mafRecords[i].Species[k-1].SLine.Start, ChromEnd: mafRecords[i].Species[k-1].SLine.Start + mafRecords[i].Species[k-1].SLine.Size, Name: "ins", Score: int64(mafRecords[i].Score)}
						bedList_del = append(bedList_del, &current_del) //append to growing bed
						bedList_ins = append(bedList_ins, &current_ins)
					}
				}
			}
		}
	}

	//write out bed files
	bed.Write(outBed_del, bedList_del, 5) //bed file has 5 fields
	bed.Write(outBed_ins, bedList_ins, 5)
	}


func usage() {
	fmt.Print(
		"humanInsertion - takes pairwise alignment maf and finds insertions in species_ins not present in species_del but flanked by continuous alignments\n" +
			"Usage:\n" +
			" humanInsertion mafFile outBed_ins outBed_del species_ins species_del\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	mafFile := flag.Arg(0)
	outBed_ins := flag.Arg(1)
	outBed_del := flag.Arg(2)
	species_ins := flag.Arg(3) //e.g. hg38
	species_del := flag.Arg(4) //e.g. panTro6

	humanInsertion(mafFile, outBed_ins, outBed_del, species_ins, species_del)
}
