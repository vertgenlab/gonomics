package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafIndels(in_maf string, species_ins string, species_del string, outIns_bed string, outDel_bed string) {
	//initialize variables
	mafRecords := maf.Read(in_maf) //Read entire in_maf. mafRecords has type Maf
	var bedList_ins []*bed.Bed     //initialize 2 bed files
	var bedList_del []*bed.Bed     //1 bed file for ins, 1 bed file for del

	//go through each line
	for i, _ := range mafRecords { //each i is a block
		for k := 1; k < len(mafRecords[i].Species); k++ { //each k is a line. Start loop at k=1, so that checking line k-1 starts at 0 and does not index out of range

			//convert maf to bed, start with getting assembly because it is needed to verify species_ins and species_del
			//here I assume only pairwise alignment, not >2 species
			//here I assume species_ins is target (because the 1st line in the block is a line, target is 2nd line in the block, k=1); species_del is query (3rd line in the block, k=2)
			assembly_del, chrom_del := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src) //start with species_del, get assembly (e.g. rheMac10), chrom (e.g. chrI)
			assembly_ins, chrom_ins := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k-1].Src) //then find corresponding species_ins line, same block, 1 line above since pairwise maf
			//TODO: to improve clarity, consider using assembly_k first, and then setting assembly_del:=assembly_k if in fact the k line is a del line

			//get eC species_del lines
			if mafRecords[i].Species[k].ELine != nil {
				if mafRecords[i].Species[k].ELine.Status == 'C' && assembly_del == species_del { //I decided to check for both Status and Src here because they are on the same data level

					//get corresponding s species_ins lines
					if assembly_ins != species_ins { //verify line k-1 is indeed species_ins
						log.Fatalf("species_ins was incorrect. Please check you have a pairwise maf file, and entered species_ins and species_del correctly") //otherwise fatal
					}
					if mafRecords[i].Species[k-1].SLine != nil { //if corresponding species_ins line is an s line

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
	bed.Write(outDel_bed, bedList_del, 5) //bed file has 5 fields
	bed.Write(outIns_bed, bedList_ins, 5)
}

func usage() {
	fmt.Print(
		"mafIndels - takes pairwise alignment maf and finds insertions in species_ins not present in species_del but flanked by continuous alignments\n" +
		"in.maf - here I assume only pairwise alignment, not >2 species\n" + //note my assumptions here
		"species_ins, species_del - here I assume species_ins is target (1st line in the block, k=0); species_del is query (2nd line in the block, k=1)\n" + //note my assumptions here
		"outIns.bed, outDel.bed - program outputs 2 bed files for species_ins and species_del coordinates, respectively. Designate filenames here\n" +
			"Usage:\n" +
			" mafIndels in.maf species_ins species_del outIns.bed outDel.bed\n" +
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

	in_maf := flag.Arg(0)
	species_ins := flag.Arg(1) //e.g. hg38
	species_del := flag.Arg(2) //e.g. rheMac10
	outIns_bed := flag.Arg(3)
	outDel_bed := flag.Arg(4)

	mafIndels(in_maf, species_ins, species_del, outIns_bed, outDel_bed)
}
