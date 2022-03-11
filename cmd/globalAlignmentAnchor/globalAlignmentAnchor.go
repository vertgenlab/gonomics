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
	"path"
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
	out_species1_filename := strings.Replace(in_maf, path.Base(in_maf), "out_"+species1+"_match.bed", 1)
	out_species2_filename := strings.Replace(in_maf, path.Base(in_maf), "out_"+species2+"_match.bed", 1)
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
				// chrom should be the same between species1 and species2 (special case: panTro6 chr2A, chr2B both count as a match to hg38 chr2)
				pT6_hg38_chr2 := ((assembly_species1 == "hg38" && assembly_species2 == "panTro6" && chrom_species1 == "chr2" && ((chrom_species2 == "chr2A") || (chrom_species2 == "chr2B"))) || (assembly_species1 == "panTro6" && assembly_species2 == "hg38" && ((chrom_species1 == "chr2A") || (chrom_species1 == "chr2B")) && chrom_species2 == "chr2")) // the special case written as a bool variable
				if chrom_species2 == chrom_species1 || pT6_hg38_chr2 {
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
	out_species1_filename := strings.Replace(in_species1_match, "match", "gap", 1)
	out_species2_filename := strings.Replace(in_species2_match, "match", "gap", 1)
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
		}

		// "else" encompasses "first entry", OR "if there is existing chr", OR the normal action "when encounter new chr"
		current_species1 = bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: species1_match_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
		current_species2 = bed.Bed{Chrom: species2_match_bed[i].Chrom, ChromStart: pos_species2, ChromEnd: species2_match_bed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}

		// before writing bed, make sure that
		// in each species, ChromStart is not equal to ChromEnd (e.g. a match entry starts at chr3 1, so the gap entry will be chr3 1 1, but can't be written to bed)
		// in each species, gap sequence should progress linearly along the chromosome (e.g. alignment match sequence skips around the chromosome, causing gap entries to skip around, ChromStart > ChromEnd)
		// the size of the gaps are practical for our alignment algorithm. The 2 sequences' product should be <=1E10. Calculate gap size product
		gapSizeProduct := (current_species1.ChromEnd - current_species1.ChromStart) * (current_species2.ChromEnd - current_species2.ChromStart)
		if !(current_species1.ChromStart < current_species1.ChromEnd && current_species2.ChromStart < current_species2.ChromEnd) {
			fmt.Printf("This bed entry pair is discarded because ChromStart or ChromEnd is invalid: %v, %v \n", current_species1, current_species2)
		} else if gapSizeProduct > 10000000000 {
			fmt.Printf("This bed entry pair is discarded because their sizes are too large: %v, %v \n", current_species1, current_species2)
		} else {
			bed.WriteBed(out_species1.File, current_species1)
			bed.WriteBed(out_species2.File, current_species2)
		}

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
	out_alignment_filename := strings.Replace(in_species1_gap, path.Base(in_species1_gap), "out_alignment.tsv", 1)
	out_alignment := fileio.EasyCreate(out_alignment_filename)
	var err error

	// initialize variables before loop
	// keep track of sequences
	var species1_seq, species2_seq []dna.Base

	// loop through input gap beds
	for i := range species1_gap_bed {
		// obtain sequences from the genome. To convert bed region (1-based, [closed,open)) to fasta index (0-based, [closed,closed]), subtract 1 from both ChromStart and ChromEnd.
		species1_seq = species1_genome_fastaMap[species1_gap_bed[i].Chrom][(species1_gap_bed[i].ChromStart - 1):(species1_gap_bed[i].ChromEnd - 1)]
		dna.AllToUpper(species1_seq) // convert all bases to uppercase, otherwise get index out of range error in scoring matrix
		species2_seq = species2_genome_fastaMap[species2_gap_bed[i].Chrom][(species2_gap_bed[i].ChromStart - 1):(species2_gap_bed[i].ChromEnd - 1)]
		dna.AllToUpper(species2_seq)

		// align with affineGap, customizeCheckersize, HumanChimpTwoScoreMatrix
		bestScore, aln := align.AffineGap_customizeCheckersize(species1_seq, species2_seq, align.HumanChimpTwoScoreMatrix, -600, -150, 10000, 10000)

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

// main function: assembles all steps
func globalAlignmentAnchor(in_maf string, species1 string, species2 string, species1_genome string, species2_genome string) {
	mafToMatch(in_maf, species1, species2)
	in_species1_match := strings.Replace(in_maf, path.Base(in_maf), "out_"+species1+"_match.bed", 1)
	in_species2_match := strings.Replace(in_maf, path.Base(in_maf), "out_"+species2+"_match.bed", 1)
	matchToGap(species1, species2, in_species1_match, in_species2_match, species1_genome, species2_genome)
	in_species1_gap := strings.Replace(in_species1_match, "match", "gap", 1)
	in_species2_gap := strings.Replace(in_species2_match, "match", "gap", 1)
	gapToAlignment(in_species1_gap, in_species2_gap, species1_genome, species2_genome)
}

func usage() {
	fmt.Print(
		"globalAlignmentAnchor - takes pairwise alignment maf, filters for trusted matches (s lines generated from the same chromosome in both species), and aligns the gap sequences between the trusted matches (affineGap, DefaultScoreMatrix)\n" +
			"in_maf - maf file. Pairwise alignment, not >2 species\n" +
			"species1, species2 - species names, e.g. hg38. Species1 is target (first line in each maf block); species2 is query (second line in each maf block)\n" +
			"species1_genome, species2_genome - fasta files containing the whole genome of each species. Each fasta sequence is 1 chromosome\n" +
			"Usage:\n" +
			"	globalAlignmentAnchor in_maf species1 species2 species1_genome species2_genome\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNum int = 5

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNum {
		flag.Usage()
		log.Fatalf("error, expecting 5 command line arguments, only found %d\n", len(flag.Args()))
	}

	in_maf := flag.Arg(0)
	species1 := flag.Arg(1)
	species2 := flag.Arg(2)
	species1_genome := flag.Arg(3)
	species2_genome := flag.Arg(4)

	globalAlignmentAnchor(in_maf, species1, species2, species1_genome, species2_genome)
}
