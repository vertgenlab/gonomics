// Command Group: "Linear Alignment Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/bed"
	"log"
	"io"
	//"os"
	"strings"
)

// helper functions to write to file
func writeToFileHandle(file io.Writer, species1 bed.Bed, species2 bed.Bed, score int64, cigar []align.Cigar) {
	var err error
	_, err = fmt.Fprintf(file, "%s\t%s\t%d\t%v\n", species1, species2, score, cigar)
	exception.PanicOnErr(err)
}

//func write(filename string, records []Bed) {
	//var err error
	//file := fileio.EasyCreate(filename)

	//for i := range records {
	//	WriteToFileHandle(file, records[i])
	//}
	//err = file.Close()
	//exception.PanicOnErr(err)
//}

// Step 1: Filter maf to remove S lines we don't trust, creating filtered maf aka anchors
// not to be confused with cmd/mafFilter, which filters for scores above a threshold
func mafToAnchor(in_maf string, species1 string, species2 string) {
	mafRecords := maf.Read(in_maf) // read input maf
	var mafFiltered []*maf.Maf // generate a container for filtered maf
	// open output bed files to write line-by-line and create variable for error
	out_species1 := fileio.EasyCreate("out_species1_match.bed")
	out_species2 := fileio.EasyCreate("out_species2_match.bed")
	var err error

	//go through each line
	//here I assume only pairwise alignment, not >2 species
	//here I assume species1 is target; species2 is query
	for i := range mafRecords { //each i is a maf block
		for k := 1; k < len(mafRecords[i].Species); k++ { //each k is a line. Start loop at k=1 because that is the lowest possible index to find species2, which is query
			//get assembly (e.g. rheMac10), chrom (e.g. chrI)
			assembly_species2, chrom_species2 := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src)
			assembly_species1, chrom_species1 := maf.SrcToAssemblyAndChrom(mafRecords[i].Species[0].Src)

			//verify line 0 is indeed species1
			if assembly_species1 != species1 {
				log.Fatalf("species1 was incorrect. Please check that you have a pairwise maf file, and entered species1 and species2 correctly")
			}

			//get s lines
			if mafRecords[i].Species[k].SLine != nil && assembly_species2 == species2 && mafRecords[i].Species[0].SLine != nil {
				//filter out only s lines that we trust to save to filtered maf
				//chrom should be the same between species1 and species2
				if chrom_species2 == chrom_species1 {
					mafFiltered = append(mafFiltered, mafRecords[i])
					// get trusted coordinates here as well. Output into bed file. Should not beyond chromosome start and end positions
					current_species2 := bed.Bed{Chrom: chrom_species2, ChromStart: mafRecords[i].Species[k].SLine.Start, ChromEnd: mafRecords[i].Species[k].SLine.Start + mafRecords[i].Species[k].SLine.Size, Name: "species2_s_filtered", Score: int(mafRecords[i].Score), FieldsInitialized: 5} //get chrom,start,end,name,score
					current_species1 := bed.Bed{Chrom: chrom_species1, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "species1_s_filtered", Score: int(mafRecords[i].Score), FieldsInitialized: 5}
					bed.WriteBed(out_species1.File, current_species1)
					bed.WriteBed(out_species2.File, current_species2)
				}
			}
		}
	}

	//generate out_maf filename and write mafFiltered to it
	out_maf := strings.Replace(in_maf, ".maf", ".filtered.maf", 1)
	maf.Write(out_maf, mafFiltered)
	// close output bed files and check for errors
	err = out_species1.Close()
	exception.PanicOnErr(err)
	err = out_species2.Close()
	exception.PanicOnErr(err)
}

//Step 2: Use anchors to calculate coordinates that still need to be aligned
func anchorToCoordinates(species1_match_bed_filename string, species2_match_bed_filename string, species1_genome_fa string, species2_genome_fa string) {
	//read bed files
	species1_match_bed := bed.Read(species1_match_bed_filename)
	species2_match_bed := bed.Read(species2_match_bed_filename)
	//read genome files, which are fastas containing each chromosome
	species1_genome_fasta := fasta.Read(species1_genome_fa)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fasta)
	species2_genome_fasta := fasta.Read(species2_genome_fa)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fasta)
	fmt.Printf("try to index species1_genome_fastaMap chr, species1_genome_fastaMap[chr1]: %v\n", species1_genome_fastaMap["chr1"])
	fmt.Printf("try to index species1_genome_fastaMap base, species1_genome_fastaMap[chr1][3]: %v\n", species1_genome_fastaMap["chr1"][3])
	fmt.Printf("try to index species1_genome_fastaMap length, len(species1_genome_fastaMap[chr1]): %v\n", len(species1_genome_fastaMap["chr1"]))
	//initialize variables to keep track of chromosome, position
	chr_prev := "" //initialize chr_prev as an empty string
	chr_curr := "" //initialize chr_curr as an empty string
	pos_species1 := 1 //initialize pos as 1. TODO: check boundaries for bed and fa. Also, maybe call from fa rather than assume start is 1
	pos_species2 := 1
	//for now, put coordinates into bed file. In the future, can just be bed object
	out_species1 := fileio.EasyCreate("out_species1_gap.bed") //rather than bedlist, write bed line by line, 1 bed for species1, 1 bed for species2
	defer out_species1.Close()
	out_species2 := fileio.EasyCreate("out_species2_gap.bed")
	defer out_species2.Close()

	//loop through bed. species1_bed and species2_match_bed should have the same number of records
	for i := range species1_match_bed {
		chr_curr = species1_match_bed[i].Chrom //set chr_curr to the new record
		//calculate the unaligned/gap chunk before we get to the aligned s line
		if i != 0 && chr_curr != chr_prev { //if this is not the first entry, but we encounter new chr
			//first finish off the previous chr
			//fmt.Printf("i, species2_match_bed[i-1].Chrom: %v, %v", i, species2_match_bed[i-1].Chrom) //TODO: remove after debugging
			current_species2 := bed.Bed{Chrom: species2_match_bed[i-1].Chrom, ChromStart: pos_species2, ChromEnd: len(species2_genome_fastaMap[species2_match_bed[i-1].Chrom]), Name: "species2_gap", FieldsInitialized: 4} //species1_genome_fa should be "FastaMap" to look up sequence name
			current_species1 := bed.Bed{Chrom: chr_prev, ChromStart: pos_species1, ChromEnd: len(species1_genome_fastaMap[chr_prev]), Name: "species1_gap", FieldsInitialized: 4}
			bed.WriteBed(out_species1.File, current_species1)
			bed.WriteBed(out_species2.File, current_species2)

			//then start the current chr
			pos_species1 = 1
			pos_species2 = 1
		} // "else" encompasses "first entry", and "if we have existing chr"
		current_species2 := bed.Bed{Chrom: species2_match_bed[i].Chrom, ChromStart: pos_species2, ChromEnd: species2_match_bed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}
		current_species1 := bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: species1_match_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
		bed.WriteBed(out_species1.File, current_species1)
		bed.WriteBed(out_species2.File, current_species2)

		//update variables at the end of each iteration
		pos_species1 = species1_match_bed[i].ChromEnd
		pos_species2 = species2_match_bed[i].ChromEnd
		chr_prev = chr_curr
		//copy is for slices, so just use equal
		//copy(pos_species1, species1_bed[i].ChromEnd)
		//copy(pos_species2, species2_match_bed[i].ChromEnd)
		//copy(chr_curr, chr_prev)

	}
	//put last entry here
	current_species2 := bed.Bed{Chrom: species2_match_bed[len(species2_match_bed)-1].Chrom, ChromStart: pos_species2, ChromEnd: len(species2_genome_fastaMap[species2_match_bed[len(species2_match_bed)-1].Chrom]), Name: "species2_gap", FieldsInitialized: 4}
	current_species1 := bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: len(species1_genome_fastaMap[chr_prev]), Name: "species1_gap", FieldsInitialized: 4}
	bed.WriteBed(out_species1.File, current_species1)
	bed.WriteBed(out_species2.File, current_species2)
}

//Step 3: globalAlignment lowMem for non-anchor sequences
//faONe is target, faTwo is query
func coordinatesToAlignment(species1_gap_bed_filename string, species2_gap_bed_filename string, species1_genome_fa string, species2_genome_fa string) {
	species1_gap_bed := bed.Read(species1_gap_bed_filename)
	species2_gap_bed := bed.Read(species2_gap_bed_filename)
	species1_genome_fasta := fasta.Read(species1_genome_fa)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fasta)
	species2_genome_fasta := fasta.Read(species2_genome_fa)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fasta)

	output := fileio.EasyCreate("out_alignment.tsv")
	var err error
	var species1_seq, species2_seq []dna.Base

	for i := range species1_gap_bed {
		species1_seq = species1_genome_fastaMap[species1_gap_bed[i].Chrom][species1_gap_bed[i].ChromStart:species1_gap_bed[i].ChromEnd] // I believe Seq starts index at 0, according to fasta.go WriteFasta function
		dna.AllToUpper(species1_seq) //convert all bases to uppercase, otherwise get index out of range error in scoring matrix
		species2_seq = species2_genome_fastaMap[species2_gap_bed[i].Chrom][species2_gap_bed[i].ChromStart:species2_gap_bed[i].ChromEnd]
		dna.AllToUpper(species2_seq)

		//alignment with affineGap
		fmt.Printf("species1_seq: %v", species1_seq)
		fmt.Printf("species2_seq: %v", species2_seq)
		bestScore, aln := align.AffineGap_customizeCheckersize(species1_seq, species2_seq, align.DefaultScoreMatrix, -400, -30, 10000, 10000)
		fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)

		//visualize
		visualize := align.View(species1_seq, species2_seq, aln)
		fmt.Println(visualize)

		//write to file
		writeToFileHandle(output, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
	}

	err = output.Close()
	exception.PanicOnErr(err)
}

//raven edited this block to specify only 1 sequnce is expected in each fasta file and add Usage nad options
func usage() {
	fmt.Print(
		"./globalAlignment - chelsea's global alignment\n" +
			" Align 2 .fasta files, each with only 1 sequence\n" +
			"Usage:\n" +
			"	globalAlignment target.fasta query.fasta\n" +
			"options:\n")
	//TODO: copied from copied from mafInspecies2s, modify
	fmt.Print(
		"mafInspecies2s - takes pairwise alignment maf and finds species1ertions in species1 not present in species2 but flanked by continuous alignments\n" +
			"in.maf - here I assume only pairwise alignment, not >2 species\n" + //note my assumptions here
			"species1, species2 - here I assume species1 is target (1st line in the block, k=0); species2 is query (2nd line in the block, k=1)\n" + //note my assumptions here
			"out_species1.bed, out_species2.bed - program outputs 2 bed files for species1 and species2 coordinates, respectively. Designate filenames here\n" +
			"Usage:\n" +
			" mafInspecies2s in.maf species1 species2 out_species1.bed out_species2.bed\n" +
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
	species1 := flag.Arg(1)
	species2 := flag.Arg(2)

	mafToAnchor(in_maf, species1, species2)
	anchorToCoordinates("testdata/species1ForStep2.bed", "testdata/species2ForStep2.bed", "testdata/species1.fa", "testdata/species2.fa")
	coordinatesToAlignment("species1ForStep3.bed", "speices2ForStep3.bed", "testdata/species1.fa", "testdata/species2.fa")
}
