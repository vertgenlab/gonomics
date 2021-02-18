package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

func mafIndels(in_maf string, species_ins string, species_del string, threshold float64, outIns_bed string, outDel_bed string) {
	//initialize variables
	mafRecords := maf.Read(in_maf) //Read entire in_maf. mafRecords has type Maf. Maf has no ReadToChan function for now
	//var bedList_ins []*bed.Bed     //initialize 2 bed files
	//var bedList_del []*bed.Bed     //1 bed file for ins, 1 bed file for del
	out_ins := fileio.EasyCreate(outIns_bed) //rather than bedlist, write bed line by line, 1 bed for ins, 1 bed for del
	defer out_ins.Close()
	out_del := fileio.EasyCreate(outDel_bed)
	defer out_del.Close() //Check if defer will work for 2 files at the same time. Seems to have worked

	//go through each line
	for i, _ := range mafRecords { //each i is a block
		for k := 1; k < len(mafRecords[i].Species); k++ { //each k is a line. Start loop at k=1 because that is the lowest possible index to find species_del, which is query

			//convert maf to bed, start with getting assembly because it is needed to verify species_ins and species_del
			//here I assume only pairwise alignment, not >2 species
			//here I assume species_ins is target (the 1st line in the block is a line but not included in mafRecords[i].Species, target is 2nd line in the block but 1st line in mafRecords[i].Species, index 0); species_del is query (3rd line in the block but 2nd line in mafRecords[i].Species, index 1 or higher)
			assembly_del, chrom_del := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src) //start with species_del, get assembly (e.g. rheMac10), chrom (e.g. chrI)
			assembly_ins, chrom_ins := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[0].Src) //then find corresponding species_ins line, which is target, index 0
			//TODO: to improve clarity, consider using assembly_k first, and then setting assembly_del:=assembly_k if in fact the k line is a del line

			//verify line 0 is indeed species_ins
			if assembly_ins != species_ins {
				log.Fatalf("species_ins was incorrect. Please check you have a pairwise maf file, and entered species_ins and species_del correctly") //otherwise fatal
			}

			//get eC species_del lines
			if mafRecords[i].Species[k].ELine != nil && assembly_del == species_del && mafRecords[i].Species[0].SLine != nil { //common checks for both eC and eI lines
				if mafRecords[i].Species[k].ELine.Status == 'C' { //Check for ELine Status here

					//convert maf to bed, continued
					current_del := bed.Bed{Chrom: chrom_del, ChromStart: mafRecords[i].Species[k].ELine.Start, ChromEnd: mafRecords[i].Species[k].ELine.Start + mafRecords[i].Species[k].ELine.Size, Name: "del_eC", Score: int64(mafRecords[i].Score)} //get chrom,start,end,name,score
					current_ins := bed.Bed{Chrom: chrom_ins, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "ins_eC", Score: int64(mafRecords[i].Score)}
					//bedList_del = append(bedList_del, &current_del) //append to growing bed
					//bedList_ins = append(bedList_ins, &current_ins)
					bed.WriteBed(out_ins.File, &current_ins, 5)
					bed.WriteBed(out_del.File, &current_del, 5)

					//get eI species_del lines
				} else if mafRecords[i].Species[k].ELine.Status == 'I' {

					//test if species_del eI fragment size < 10% corresponding s fragment size
					//make sure arithmetic is all on float64
					if float64(mafRecords[i].Species[k].ELine.Size) < threshold*float64(mafRecords[i].Species[0].SLine.Size) {

						//convert maf to bed, continued
						current_del := bed.Bed{Chrom: chrom_del, ChromStart: mafRecords[i].Species[k].ELine.Start, ChromEnd: mafRecords[i].Species[k].ELine.Start + mafRecords[i].Species[k].ELine.Size, Name: "del_eI", Score: int64(mafRecords[i].Score)} //get chrom,start,end,name,score
						current_ins := bed.Bed{Chrom: chrom_ins, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "ins_eI", Score: int64(mafRecords[i].Score)}
						//bedList_del = append(bedList_del, &current_del) //append to growing bed
						//bedList_ins = append(bedList_ins, &current_ins)
						bed.WriteBed(out_ins.File, &current_ins, 5)
						bed.WriteBed(out_del.File, &current_del, 5)
					}
				}
			}
		}
	}
	//write out bed files
	//bed.Write(outDel_bed, bedList_del, 5) //bed file has 5 fields
	//bed.Write(outIns_bed, bedList_ins, 5)
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
	eiThreshold := flag.Float64("eiThreshold", 0.1, "fraction of insertion size that marks the maximum acceptable threshold of unaligned fragment in species_del")
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
	threshold := *eiThreshold

	mafIndels(in_maf, species_ins, species_del, threshold, outIns_bed, outDel_bed)
}
