// Command Group: "Linear Alignment Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
	"io"
	"log"
	"strings"
)

// helper function: write to tsv file
func writeToFileHandle(file io.Writer, species1 bed.Bed, species2 bed.Bed, score int64, cigar []align.Cigar) {
	var err error
	_, err = fmt.Fprintf(file, "%s\t%s\t%d\t%v\n", species1, species2, score, cigar)
	exception.PanicOnErr(err)
}

// Step 1: Filter maf to remove S lines we don't trust, creating filtered maf (aka anchors, or "match")
// not to be confused with cmd/mafFilter, which filters for scores above a threshold
func mafToMatch(in_maf string, species1 string, species2 string) {
	mafRecords := maf.Read(in_maf) // read input maf

	// open output files to write line-by-line and create variable for error
	out_maf_filename := strings.Replace(in_maf, ".maf", ".filtered.maf", 1)
	out_species1_filename := strings.Replace("out_species1_match.bed", "species1", species1, 1)
	out_species2_filename := strings.Replace("out_species2_match.bed", "species2", species2, 1)
	out_maf := fileio.EasyCreate(out_maf_filename)
	out_species1 := fileio.EasyCreate(out_species1_filename)
	out_species2 := fileio.EasyCreate(out_species2_filename)
	var err error

	// initialize variables before loop
	// keep track of assembly, chromosome
	var assembly_species1, assembly_species2, chrom_species1, chrom_species2 string
	// containers for entries to write to ouput files
	var bed_species1, bed_species2 bed.Bed

	// loop through input maf
	// I assume pairwise alignment, not >2 species
	// I assume species1 is target; species2 is query
	for i := range mafRecords { // each i is a maf block
		// get assembly (e.g. rheMac10), chrom (e.g. chrI)
		assembly_species1, chrom_species1 = maf.SrcToAssemblyAndChrom(mafRecords[i].Species[0].Src)
		// get trusted match coordinates here as well. Output into bed file. Get chrom,start,end,name,score
		bed_species1 = bed.Bed{Chrom: chrom_species1, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "species1_s_filtered_match", Score: int(mafRecords[i].Score), FieldsInitialized: 5}

		for k := 1; k < len(mafRecords[i].Species); k++ { // each k is a line. Start loop at k=1 because that is the lowest possible index to find species2, which is query
			assembly_species2, chrom_species2 = maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src)

			// verify line 0 is indeed species1
			if assembly_species1 != species1 {
				log.Fatalf("species1 was incorrect. Please check that you have a pairwise maf file, and entered species1 and species2 correctly")
			}

			// get s lines
			if mafRecords[i].Species[k].SLine != nil && assembly_species2 == species2 && mafRecords[i].Species[0].SLine != nil {

				// filter out only s lines that we trust to save to filtered maf
				// chrom should be the same between species1 and species2
				if chrom_species2 == chrom_species1 {
					maf.WriteToFileHandle(out_maf, mafRecords[i])
					bed_species2 = bed.Bed{Chrom: chrom_species2, ChromStart: mafRecords[i].Species[k].SLine.Start, ChromEnd: mafRecords[i].Species[k].SLine.Start + mafRecords[i].Species[k].SLine.Size, Name: "species2_s_filtered_match", Score: int(mafRecords[i].Score), FieldsInitialized: 5}
					bed.WriteBed(out_species1.File, bed_species1)
					bed.WriteBed(out_species2.File, bed_species2)
				}
			}
		}
	}

	// close output files and check for errors
	err = out_maf.Close()
	exception.PanicOnErr(err)
	err = out_species1.Close()
	exception.PanicOnErr(err)
	err = out_species2.Close()
	exception.PanicOnErr(err)
}

// Step 2: Use match to calculate coordinates that still need to be aligned, aka "gap"
func matchToGap(species1 string, species2 string, in_species1_match string, in_species2_match string, species1_genome string, species2_genome string) {
	// read input files
	species1_match_bed := bed.Read(in_species1_match)
	species2_match_bed := bed.Read(in_species2_match)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_species1_filename := strings.Replace("out_species1_gap.bed", "species1", species1, 1)
	out_species2_filename := strings.Replace("out_species2_gap.bed", "species2", species2, 1)
	out_species1 := fileio.EasyCreate(out_species1_filename)
	out_species2 := fileio.EasyCreate(out_species2_filename)
	var err error

	// initialize variables before loop
	// keep track of chromosome, position
	chr_prev := "" // initialize chr as an empty string
	chr_curr := ""
	pos_species1 := 1 // initialize pos as 1. bed and fa both start at 1
	pos_species2 := 1
	// containers for entries to write to ouput files
	var current_species1, current_species2 bed.Bed

	// loop through input match beds
	// species1_match_bed and species2_match_bed should have the same number of records, and can be indexed simultaneously
	for i := range species1_match_bed {
		chr_curr = species1_match_bed[i].Chrom // set chr_curr to the new record

		// calculate the unaligned/gap region before arriving at the aligned/match s line
		if i != 0 && chr_curr != chr_prev { // if this is not the first entry, AND encounter new chr, additional action is required
			// first finish off the previous chr
			current_species1 = bed.Bed{Chrom: chr_prev, ChromStart: pos_species1, ChromEnd: len(species1_genome_fastaMap[chr_prev]), Name: "species1_gap", FieldsInitialized: 4}
			current_species2 = bed.Bed{Chrom: species2_match_bed[i-1].Chrom, ChromStart: pos_species2, ChromEnd: len(species2_genome_fastaMap[species2_match_bed[i-1].Chrom]), Name: "species2_gap", FieldsInitialized: 4}
			bed.WriteBed(out_species1.File, current_species1)
			bed.WriteBed(out_species2.File, current_species2)

			// then start the current chr
			pos_species1 = 1
			pos_species2 = 1
		} // "else" encompasses "first entry", OR "if there is existing chr", OR the normal action "when encounter new chr"
		current_species1 = bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: species1_match_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
		current_species2 = bed.Bed{Chrom: species2_match_bed[i].Chrom, ChromStart: pos_species2, ChromEnd: species2_match_bed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}
		bed.WriteBed(out_species1.File, current_species1)
		bed.WriteBed(out_species2.File, current_species2)

		// update variables at the end of each iteration
		pos_species1 = species1_match_bed[i].ChromEnd
		pos_species2 = species2_match_bed[i].ChromEnd
		chr_prev = chr_curr
	}

	// after loop, need to write the last entry to output files
	current_species1 = bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: len(species1_genome_fastaMap[chr_prev]), Name: "species1_gap", FieldsInitialized: 4}
	current_species2 = bed.Bed{Chrom: species2_match_bed[len(species2_match_bed)-1].Chrom, ChromStart: pos_species2, ChromEnd: len(species2_genome_fastaMap[species2_match_bed[len(species2_match_bed)-1].Chrom]), Name: "species2_gap", FieldsInitialized: 4}
	bed.WriteBed(out_species1.File, current_species1)
	bed.WriteBed(out_species2.File, current_species2)

	// close output files and check for errors
	err = out_species1.Close()
	exception.PanicOnErr(err)
	err = out_species2.Close()
	exception.PanicOnErr(err)
}

// Step 3: align "gap" sequences
func gapToAlignment(in_species1_gap string, in_species2_gap string, species1_genome string, species2_genome string) {
	// read input files
	species1_gap_bed := bed.Read(in_species1_gap)
	species2_gap_bed := bed.Read(in_species2_gap)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_alignment := fileio.EasyCreate("out_alignment.tsv")
	var err error

	// initialize variables before loop
	// keep track of sequences
	var species1_seq, species2_seq []dna.Base

	// loop through input gap beds
	for i := range species1_gap_bed {
		// obtain sequences from the genome
		species1_seq = species1_genome_fastaMap[species1_gap_bed[i].Chrom][species1_gap_bed[i].ChromStart:species1_gap_bed[i].ChromEnd]
		dna.AllToUpper(species1_seq) // convert all bases to uppercase, otherwise get index out of range error in scoring matrix
		species2_seq = species2_genome_fastaMap[species2_gap_bed[i].Chrom][species2_gap_bed[i].ChromStart:species2_gap_bed[i].ChromEnd]
		dna.AllToUpper(species2_seq)

		// align with affineGap, customizeCheckersize, DefaultScoreMatrix
		bestScore, aln := align.AffineGap_customizeCheckersize(species1_seq, species2_seq, align.DefaultScoreMatrix, -400, -30, 10000, 10000)

		// optional: print results to terminal
		// print alignment score, cigar
		//fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)
		// visualize
		//visualize := align.View(species1_seq, species2_seq, aln)
		//fmt.Println(visualize)

		// write to output file
		writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
	}

	// close output files and check for errors
	err = out_alignment.Close()
	exception.PanicOnErr(err)
}

//TODO: start editing here
//raven edited this block to specify only 1 sequnce is expected in each fasta file and add Usage nad options
func usage() {
	fmt.Print(
		"./globalAlignment - chelsea's global alignment\n" +
			" Align 2 .fasta files, each with only 1 sequence\n" +
			"Usage:\n" +
			"	globalAlignment target.fasta query.fasta\n" +
			"options:\n")
	//TODO: copied from copied from mafInspecies2s, modify
	//note whole genome fasta with all chromosomes
	// align "gap" sequences with affineGap, customizeCheckersize, DefaultScoreMatrix
	fmt.Print(
		"mafInspecies2s - takes pairwise alignment maf and finds species1ertions in species1 not present in species2 but flanked by continuous alignments\n" +
			"in.maf - here I assume only pairwise alignment, not >2 species\n" +
			"species1, species2 - here I assume species1 is target (1st line in the block, k=0); species2 is query (2nd line in the block, k=1)\n" +
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

	mafToMatch(in_maf, species1, species2)
	matchToGap(species1, species2, "testdata/out_hg38_match.bed", "testdata/out_rheMac10_match.bed", "testdata/species1.fa", "testdata/species2.fa")
	gapToAlignment("testdata/out_hg38_gap.bed", "testdata/out_rheMac10_gap.bed", "testdata/species1.fa", "testdata/species2.fa")
}
