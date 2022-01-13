// Command Group: "Linear Alignment Tools"

package main

import (
	"flag"
	"fmt"
	//"github.com/vertgenlab/gonomics/align"
	//"github.com/vertgenlab/gonomics/exception"
	//"github.com/vertgenlab/gonomics/fasta"
	//"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/maf"
	"log"
	//"os"
	"strings"
)

//Step 1: Filter maf to remove S lines we don't trust, creating filtered maf aka anchors
//not to be confused with cmd/mafFilter, which filters for scores above a threshold
func mafToAnchor(in_maf string, species_ins string, species_del string) {
	//read in_maf and generate a container for mafFiltered
	mafRecords := maf.Read(in_maf)
	var mafFiltered []*maf.Maf

	//go through each line
	//here I assume only pairwise alignment, not >2 species
	//here I assume species_ins is target; species_del is query
	for i := range mafRecords { //each i is a maf block
		for k := 1; k < len(mafRecords[i].Species); k++ { //each k is a line. Start loop at k=1 because that is the lowest possible index to find species_del, which is query
			//get assembly (e.g. rheMac10), chrom (e.g. chrI)
			assembly_del, chrom_del := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src)
			assembly_ins, chrom_ins := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[0].Src)

			//verify line 0 is indeed species_ins
			if assembly_ins != species_ins {
				log.Fatalf("species_ins was incorrect. Please check that you have a pairwise maf file, and entered species_ins and species_del correctly")
			}

			//get s lines
			if mafRecords[i].Species[k].SLine != nil && assembly_del == species_del && mafRecords[i].Species[0].SLine != nil {
				//filter out only s lines that we trust to save to filtered maf
				//chrom should be the same between species_ins and species_del
				if chrom_del == chrom_ins {
					mafFiltered = append(mafFiltered, mafRecords[i])
				}
			}
		}
	}

	//generate out_maf filename and write mafFiltered to it
	out_maf := strings.Replace(in_maf, ".maf", ".filtered.maf", 1)
	maf.Write(out_maf, mafFiltered)
}

//Step 2: Use anchors aka filtered maf to calculate coordinates that still need to be aligned
//TODO: save as bed entry, but what to do about strand
func anchorToCoordinates(in_maf string, species_ins string, species_del string, ins_genome_fa string, del_genome_fa string) {
	//initialize variables
	mafRecords := maf.Read(in_maf) //Read entire in_maf. mafRecords has type Maf. Maf has no ReadToChan function for now
	//var bedList_ins []*bed.Bed     //initialize 2 bed files
	//var bedList_del []*bed.Bed     //1 bed file for ins, 1 bed file for del

	//started work for anchorToCoordinates here
	//template is cmd/faFilter
	//read genome files, which are fastas
	ins_genome := fasta.Read(ins_genome_fa)
	del_genome := fasta.Read(del_genome_fa)
	for i := 0; i < len(ins_genome); i++ {
		if records[i].Name == name { //replace "name" with chr
			return records[i].Seq[start:end] //replace "start:end" with positions that are not s. I believe Seq starts index at 0, according to fasta.go WriteFasta function
		}
	}

	out_ins := fileio.EasyCreate(outIns_bed) //rather than bedlist, write bed line by line, 1 bed for ins, 1 bed for del
	defer out_ins.Close()
	out_del := fileio.EasyCreate(outDel_bed)
	defer out_del.Close() //Check if defer will work for 2 files at the same time. Seems to have worked

	//go through each line
	for i := range mafRecords { //each i is a block
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
					current_del := bed.Bed{Chrom: chrom_del, ChromStart: mafRecords[i].Species[k].ELine.Start, ChromEnd: mafRecords[i].Species[k].ELine.Start + mafRecords[i].Species[k].ELine.Size, Name: "del_eC", Score: int(mafRecords[i].Score), FieldsInitialized: 5} //get chrom,start,end,name,score
					current_ins := bed.Bed{Chrom: chrom_ins, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "ins_eC", Score: int(mafRecords[i].Score), FieldsInitialized: 5}
					//bedList_del = append(bedList_del, &current_del) //append to growing bed
					//bedList_ins = append(bedList_ins, &current_ins)
					bed.WriteBed(out_ins.File, current_ins)
					bed.WriteBed(out_del.File, current_del)

					//get eI species_del lines
				} else if mafRecords[i].Species[k].ELine.Status == 'I' {

					//test if species_del eI fragment size < 10% corresponding s fragment size
					//make sure arithmetic is all on float64
					if float64(mafRecords[i].Species[k].ELine.Size) < threshold*float64(mafRecords[i].Species[0].SLine.Size) {

						//convert maf to bed, continued
						current_del := bed.Bed{Chrom: chrom_del, ChromStart: mafRecords[i].Species[k].ELine.Start, ChromEnd: mafRecords[i].Species[k].ELine.Start + mafRecords[i].Species[k].ELine.Size, Name: "del_eI", Score: int(mafRecords[i].Score), FieldsInitialized: 5} //get chrom,start,end,name,score
						current_ins := bed.Bed{Chrom: chrom_ins, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "ins_eI", Score: int(mafRecords[i].Score), FieldsInitialized: 5}
						//bedList_del = append(bedList_del, &current_del) //append to growing bed
						//bedList_ins = append(bedList_ins, &current_ins)
						bed.WriteBed(out_ins.File, current_ins)
						bed.WriteBed(out_del.File, current_del)
					}
				}
			}
		}
	}
	//write out bed files
	//bed.Write(outDel_bed, bedList_del, 5) //bed file has 5 fields
	//bed.Write(outIns_bed, bedList_ins, 5)
}

//Step 3: globalAlignment lowMem for non-anchor sequences
//raven did not put this helper function into the globalAlignment function because it is used twice within the globalAlignment function
//raven wrote this block to count sequences based on the Read function in gonomics/fasta/fasta.go
//raven changed the input variable from filename string to inputFile EasyReader, so that the file is only opened 1 time for 2 purposes: faDone and CountSeqIdx
func CountSeqIdx(inputFile *fileio.EasyReader) int {
	var line string
	var seqIdx int = 1 //I know in Read seqIdx is int64 and starts with -1, but I am using it differently here. EasyReader comes in having read the first fasta, so seqIdx starts with 1
	var doneReading bool = false
	for line, doneReading = fileio.EasyNextRealLine(inputFile); !doneReading; line, doneReading = fileio.EasyNextRealLine(inputFile) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
		}
	}
	err := inputFile.Close()
	exception.PanicOnErr(err)
	return seqIdx
}

//raven did not put this helper function into the globalAlignment function because it is used in the globalAlignment function
//faONe is target, faTwo is query
func cigarToGraph(target fasta.Fasta, query fasta.Fasta, aln []align.Cigar) *genomeGraph.GenomeGraph {
	answer := genomeGraph.EmptyGraph()
	//use targetEnd and queryEnd to track position number for each fasta sequence as we move along.
	var targetEnd, queryEnd int = 0, 0
	var curr *genomeGraph.Node
	//TODO: Make sure it can handle non-match at node 0. Fill out annotations for each node (knowing which allele/haplotype).
	//TODO: Handle multi-entry fasta?

	//Creating the first node. This is done independent of all other number because this node has no 'previous' nodes.  All following code leverages the cigar output (second number printed in the Op position of the struct {}) to determine if the alignment returned a 0:match, 1:insertion, or 2:deletion. All indels are relative to the target sequence.
	curr, targetEnd, queryEnd = genomeGraph.FaSeqToNode(target, query, targetEnd, queryEnd, aln[0], 0)
	curr = genomeGraph.AddNode(answer, curr)
	//Drawing the remaining nodes and all edges. Method for adding edges is based on previous nodes.
	for i := 1; i < len(aln); i++ {
		curr, targetEnd, queryEnd = genomeGraph.FaSeqToNode(target, query, targetEnd, queryEnd, aln[i], i)
		curr = genomeGraph.AddNode(answer, curr)
		if aln[i].Op == align.ColM {
			genomeGraph.AddEdge(&answer.Nodes[i-1], curr, 1)
			genomeGraph.AddEdge(&answer.Nodes[i-2], curr, 0.5)
		} else if aln[i].Op == align.ColI || aln[i].Op == align.ColD {
			genomeGraph.AddEdge(&answer.Nodes[i-1], curr, 0.5)
		} else {
			log.Fatalf("Error: cigar.Op = %d is unrecognized...\n", aln[i].Op)
		}
	}
	return answer
}

//raven moved helper functions from main and non-usage functions into this function
func globalAlignment(inputFileOne *fileio.EasyReader, inputFileTwo *fileio.EasyReader, outFileName string) {
	//make sure files meet the usage requirements of globalAlignment.go
	faOne, faDoneOne := fasta.NextFasta(inputFileOne)
	faTwo, faDoneTwo := fasta.NextFasta(inputFileTwo)
	if faDoneOne || faDoneTwo {
		log.Fatalf("error, unable to read .fa files, check for > symbol before each name and that each fasta entry is only one line. Check to make sure there are only four bases or N, globalAlignment is unable to use anything outside A-T-G-C-N.\n")
	}
	numSeqOne := CountSeqIdx(inputFileOne)
	numSeqTwo := CountSeqIdx(inputFileTwo)
	if numSeqOne > 1 || numSeqTwo > 1 {
		log.Fatalf("multiple sequnces detected in .fa files: %v sequences in the first .fa file and %v sequences in the second .fa file. This program is designed for .fa files with only 1 sequence in them\n", numSeqOne, numSeqTwo)
	}

	//needleman wunsch (global alignment)
	bestScore, aln := align.ConstGap(faOne.Seq, faTwo.Seq, align.HumanChimpTwoScoreMatrix, -430)
	fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)

	//visualize
	visualize := align.View(faOne.Seq, faTwo.Seq, aln)
	fmt.Println(visualize)

	//raven added this block to put visualized alignment output into MSA Fasta
	if outFileName != "" { //only write to file if a filename is given in the faOut command line option
		outFile, err := os.Create(outFileName) //fileio.EasyCreate and other fileio tools should work, but I don't really know how to use them
		if err != nil {
			log.Fatalf("Write to file failed on step 1\n")
		}
		visualizeOutput := ">" + faOne.Name + "\n" + strings.Split(visualize, "\n")[0] + "\n" + ">" + faTwo.Name + "\n" + strings.Split(visualize, "\n")[1] + "\n"
		_, err = outFile.WriteString(visualizeOutput)
		if err != nil {
			log.Fatalf("Write to file failed on step 2\n")
		}
		outFile.Close() //commented out defer outFile.Close()
	}

	//cigar to graph
	//genomeGraph := cigarToGraph(faOne, faTwo, aln)
	//genomeGraph.PrintGraph(genomeGraph)
}


//raven edited this block to specify only 1 sequnce is expected in each fasta file and add Usage nad options
func usage() {
	fmt.Print(
		"./globalAlignment - chelsea's global alignment\n" +
			" Align 2 .fasta files, each with only 1 sequence\n" +
			"Usage:\n" +
			"	globalAlignment target.fasta query.fasta\n" +
			"options:\n")
	//TODO: copied from copied from mafIndels, modify
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
	var expectedNum int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	//faOut := flag.String("faOut", "", "fasta MSA output filename")
	flag.Parse()

	if len(flag.Args()) != expectedNum {
		flag.Usage()
		log.Fatalf("error, expecting 2 .fasta files to be able to align, only found %d files...\n", len(flag.Args()))
	}

	//read in sequences that should be put in as fasta type files.
	//raven edited this block to save fileio.EasyOpen as file handles, so that the file is only opened 1 time for 2 purposes: faDone and CountSeqIdx
	in_maf := flag.Arg(0)
	species_ins := flag.Arg(1)
	species_del := flag.Arg(2)

	mafToAnchor(in_maf, species_ins, species_del)
}
