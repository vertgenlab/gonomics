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

// helper function: make chr name match map
func makeChrMap(chrMap_filename string) map[string][]string {

	chrMap := make(map[string][]string) // map to hold species1 species2 chr name matches
	var chrMap_stringSplit []string
	var exists bool

	chrMap_string := fileio.Read(chrMap_filename) // read file so each line is 1 string

	for i := range chrMap_string {
		chrMap_stringSplit = strings.Split(chrMap_string[i], "\t") // convert each line = 1 string into a slice of 2 strings (assumption: each line only has 2 columns separated by 1 tab). Note that Notepad generates proper tab, but not Atom
		_, exists = chrMap[chrMap_stringSplit[0]]                  // convert each line = 1 slice of 2 strings into map (bypass struct), so key is species1 chr name = slice[0] and value is a slice with all matching species2 chr name = slice[1] in >=1 lines
		if !exists {
			chrMap_species2 := make([]string, 1)       // slice to hold just 1 species2 chr name the first time its species1 chr name appears
			chrMap_species2[0] = chrMap_stringSplit[1] // use species2 chr name to populate slice of length 1 in preparation for [string][]string key-value pair. This slice will be overwritten each time there is a new species2 chr name
			//TODO: Can I define chrMap_species2 []string outside of for loop? The below line doesn't work. Because of this, if I define chrMap_species2 []string outside of for loop, I can't use "copy" to prevent the chrMap map from always pointing to the chrMap_species2 []string when chrMap_species2 gets updated
			//copy(chrMap[chrMap_stringSplit[0]], chrMap_species2)
			chrMap[chrMap_stringSplit[0]] = chrMap_species2
		} else {
			chrMap[chrMap_stringSplit[0]] = append(chrMap[chrMap_stringSplit[0]], chrMap_stringSplit[1]) // if species1 chr name already exists, append to the exisitng species2 chr name slice the new species2 chr name
		}
	}

	return chrMap
}

// helper function: check if maf entry passes checks
func matchMafPass(assembly_species1 string, assembly_species2 string, chrom_species1 string, chrom_species2 string, species1_SrcSize int, species2_SrcSize int, species1_ChromStart int, species1_ChromEnd int, species2_ChromStart int, species2_ChromEnd int, chrMap map[string][]string, diagonal bool) bool {

	pass := true

	// chrom should match between species1 and species2
	// for each chrom_species2, go through chrMap to see if it's contained within chrom_species1's matching species2 chr names
	chrMatch := false // separate bool is needed for chrMatch, start out as false
	for _, s := range chrMap[chrom_species1] {
		if chrom_species2 == s {
			chrMatch = true // only become true if chrMatch
		}
	}
	if !chrMatch { // if chrMatch == false, then pass must be set to false
		pass = false
	}

	// maf entry should be roughly diagonal
	// unless user specifies that they do not want the diagonal check
	if diagonal {
		if (float64(species2_ChromStart) <= float64(species1_ChromStart)-0.05*float64(species1_SrcSize)) || (float64(species2_ChromStart) >= float64(species1_ChromStart)+0.05*float64(species1_SrcSize)) {
			pass = false
		} else if (float64(species1_ChromStart) <= float64(species2_ChromStart)-0.05*float64(species2_SrcSize)) || (float64(species1_ChromStart) >= float64(species2_ChromStart)+0.05*float64(species2_SrcSize)) {
			pass = false
		}
	}

	return pass
}

// helper function: check if gap bed entry passes checks
// need both pos_species (most recently updated alignment position)
// and species1_ChromStart (the start of the current bed entry)
func gapBedPass(pos_species1 int, species1_ChromStart int, species1_ChromEnd int, pos_species2 int, species2_ChromStart int, species2_ChromEnd int, gapSizeProductLimit int) (bool, string, string) {
	pass := true
	species1_Name := "species1_gap"
	species2_Name := "species2_gap"
	// calculate gap size of current bed entry, based on species1_ChromStart
	species1_gapSize := species1_ChromEnd - species1_ChromStart
	species2_gapSize := species2_ChromEnd - species2_ChromStart
	// calculate gap size since the most recently updated alignment position, based on pos_species
	species1_gapSizeBig := species1_ChromEnd - pos_species1
	species2_gapSizeBig := species2_ChromEnd - pos_species2
	// define limits and necessary calculations to evaluate limit
	// gap size product limit should be evaluted based on gapSize
	gapSizeProduct := species1_gapSize * species2_gapSize
	// gap size multiple should be evaluated based on gapSizeBig
	gapSizeBigMultiple := 0.0
	if species1_gapSizeBig != 0 {
		gapSizeBigMultiple = float64(species2_gapSizeBig / species1_gapSizeBig)
	}
	gapSizeBigMultipleLimit := 100.00 // gapSizeBigMultipleLimit is currently hardcoded tentatively, 100

	if species1_gapSize > 0 && species2_gapSize == 0 { // insertion in species1, no need for alignment
		species1_Name = "species1_Insertion"
		species2_Name = "species2_gap_size0"
	} else if species1_gapSize == 0 && species2_gapSize > 0 { // insertion in species2, no need for alignment
		species2_Name = "species2_Insertion"
		species1_Name = "species1_gap_size0"
	} else if !(species1_gapSize > 0 && species2_gapSize > 0) { // need to check, otherwise may have chromEnd<chromStart, leading to gapSize<0. This is one way to make sure in each species esp the non-reference, gap sequence should progress linearly along the chromosome (e.g. alignment match sequence skips around the chromosome, causing gap entries to skip around, ChromStart > ChromEnd)
		pass = false
		species1_Name = "species1_gap,doNotCalculate_invalidChromStartOrChromEnd"
		species2_Name = "species2_gap,doNotCalculate_invalidChromStartOrChromEnd"
	} else if gapSizeBigMultiple > gapSizeBigMultipleLimit { // This is one way to make sure in each species esp the non-reference, gap sequence should progress linearly along the chromosome
		// If the currently aligned sequence in species2 is really far away from the previously aligned sequence in species2,
		// 100 times further than the currently aligned sequence in species1 is from the previously aligned sequence in species1,
		// then the current alignment should be discarded and not "trusted" to be a true alignment,
		// because the currently aligned sequence in species2 is probably a random region further downstream on the chromosome, which breaks synteny.
		// This is only a concern for species2 because species1 is always progressing linearly (has synteny) and a big gap in species1 might just mean a big insertion in species1.
		pass = false
		species1_Name = "species1_gap,doNotCalculate_largeGapSizeMultiple"
		species2_Name = "species2_gap,doNotCalculate_largeGapSizeMultiple"
		// still check for diagonal. If diagonal, can accept (rescue), but add label
		// this step should come last, after initially accepting some gaps (without aligning), rejecting other gaps, and finally in this step rescuing some rejected gaps to accept them again
		if float64(species2_ChromStart) >= 0.95*float64(species1_ChromStart) && float64(species2_ChromStart) <= 1.05*float64(species1_ChromStart) && float64(species2_ChromEnd) >= 0.95*float64(species1_ChromEnd) && float64(species2_ChromEnd) <= 1.05*float64(species2_ChromEnd) {
			pass = true
			species1_Name = "species1_gap_largeGapSize_diagonal"
			species2_Name = "species2_gap_largeGapSize_diagonal"
		}
	}

	if gapSizeProduct > gapSizeProductLimit { // the product of the gap sizes needs to be practical for our alignment algorithm. The 2 sequences' product should be <=1E10
		pass = false
		species1_Name = species1_Name + ",doNotCalculate_largeGapSizeProduct"
		species2_Name = species2_Name + ",doNotCalculate_largeGapSizeProduct"
	}

	return pass, species1_Name, species2_Name
}

// Step 1: Filter maf to remove S lines we don't trust, creating filtered maf (aka anchors, or "match")
// not to be confused with cmd/mafFilter, which filters for scores above a threshold
func mafToMatch(in_maf string, species1 string, species2 string, out_filename_prefix string, chrMap_filename string, diagonal bool) (string, string) {
	mafRecords := maf.Read(in_maf) // read input maf
	chrMap := makeChrMap(chrMap_filename)

	// open output files to write line-by-line and create variable for error
	out_maf_filename := out_filename_prefix + ".filtered.maf"
	out_species1_filename := out_filename_prefix + "_" + species1 + "_match.bed"
	out_species2_filename := out_filename_prefix + "_" + species2 + "_match.bed"
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

				bed_species2 = bed.Bed{Chrom: chrom_species2, ChromStart: mafRecords[i].Species[k].SLine.Start, ChromEnd: mafRecords[i].Species[k].SLine.Start + mafRecords[i].Species[k].SLine.Size, Name: "species2_s_filtered_match", Score: int(mafRecords[i].Score), FieldsInitialized: 5}

				// filter out only s lines that we trust to save to filtered maf
				pass := matchMafPass(assembly_species1, assembly_species2, chrom_species1, chrom_species2, mafRecords[i].Species[0].SLine.SrcSize, mafRecords[i].Species[k].SLine.SrcSize, bed_species1.ChromStart, bed_species1.ChromEnd, bed_species2.ChromStart, bed_species2.ChromEnd, chrMap, diagonal)
				if pass {
					maf.WriteToFileHandle(out_maf, mafRecords[i])
					bed.WriteBed(out_species1, bed_species1)
					bed.WriteBed(out_species2, bed_species2)
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

	// return output filenames
	return out_species1_filename, out_species2_filename
}

// Step 2: Use match to calculate coordinates that still need to be aligned, aka "gap"
func matchToGap(in_species1_match string, in_species2_match string, species1_genome string, species2_genome string, species1 string, species2 string, gapSizeProductLimit int, out_filename_prefix string) (string, string) {
	// read input files
	species1_match_bed := bed.Read(in_species1_match)
	species2_match_bed := bed.Read(in_species2_match)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_species1_filename := out_filename_prefix + "_" + species1 + "_gap.bed"
	out_species2_filename := out_filename_prefix + "_" + species2 + "_gap.bed"
	out_species1_doNotCalculate_filename := out_filename_prefix + "_" + species1 + "_gap_doNotCalculate.bed"
	out_species2_doNotCalculate_filename := out_filename_prefix + "_" + species2 + "_gap_doNotCalculate.bed"
	out_species1 := fileio.EasyCreate(out_species1_filename)
	out_species2 := fileio.EasyCreate(out_species2_filename)
	out_species1_doNotCalculate := fileio.EasyCreate(out_species1_doNotCalculate_filename)
	out_species2_doNotCalculate := fileio.EasyCreate(out_species2_doNotCalculate_filename)
	var err error

	// initialize variables before loop
	// keep track of chromosome, position
	chr_prev_species1 := species1_match_bed[0].Chrom // initialize prev chr to be the same as curr chr, so the first entry won't have "chr_curr_species1 != chr_prev_species1_species1" and be treated like changing chr
	chr_curr_species1 := species1_match_bed[0].Chrom // initialize curr chr
	chr_prev_species2 := species2_match_bed[0].Chrom
	chr_curr_species2 := species2_match_bed[0].Chrom
	pos_species1 := 1 // initialize pos as 1. bed and fa both start at 1
	pos_species2 := 1
	// check for gapBedPass
	var pass bool
	// containers for entries to write to ouput files
	var current_species1, current_species2 bed.Bed

	// write i==0 case outside of loop to avoid checking for i==0 in each iteration
	current_species1 = bed.Bed{Chrom: chr_curr_species1, ChromStart: pos_species1, ChromEnd: species1_match_bed[0].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
	current_species2 = bed.Bed{Chrom: chr_curr_species2, ChromStart: pos_species2, ChromEnd: species2_match_bed[0].ChromStart, Name: "species2_gap", FieldsInitialized: 4}
	// before writing bed, test if the bed entries pass filtering criteria
	pass, current_species1.Name, current_species2.Name = gapBedPass(pos_species1, current_species1.ChromStart, current_species1.ChromEnd, pos_species2, current_species2.ChromStart, current_species2.ChromEnd, gapSizeProductLimit)
	if !pass {
		bed.WriteBed(out_species1_doNotCalculate, current_species1)
		bed.WriteBed(out_species2_doNotCalculate, current_species2)
	} else {
		bed.WriteBed(out_species1, current_species1)
		bed.WriteBed(out_species2, current_species2)
		// update variables at the end of each iteration
		// only update position if bed passes
		pos_species1 = species1_match_bed[0].ChromEnd
		pos_species2 = species2_match_bed[0].ChromEnd
	}

	// starting from i==1
	// loop through input match beds
	// species1_match_bed and species2_match_bed should have the same number of records, and can be indexed simultaneously
	for i := 1; i < len(species1_match_bed); i++ {
		chr_curr_species1 = species1_match_bed[i].Chrom // set chr_curr_species1 to the new record
		chr_curr_species2 = species2_match_bed[i].Chrom

		// calculate the unaligned/gap region before arriving at the aligned/match s line
		if chr_curr_species1 != chr_prev_species1 { // if encounter new chr
			// first finish off the previous chr
			current_species1 = bed.Bed{Chrom: chr_prev_species1, ChromStart: species1_match_bed[i-1].ChromEnd, ChromEnd: len(species1_genome_fastaMap[chr_prev_species1]), Name: "species1_gap", FieldsInitialized: 4}
			current_species2 = bed.Bed{Chrom: chr_prev_species2, ChromStart: species2_match_bed[i-1].ChromEnd, ChromEnd: len(species2_genome_fastaMap[chr_prev_species2]), Name: "species2_gap", FieldsInitialized: 4}
			pass, current_species1.Name, current_species2.Name = gapBedPass(pos_species1, current_species1.ChromStart, current_species1.ChromEnd, pos_species2, current_species2.ChromStart, current_species2.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(out_species1_doNotCalculate, current_species1)
				bed.WriteBed(out_species2_doNotCalculate, current_species2)
			} else {
				bed.WriteBed(out_species1, current_species1)
				bed.WriteBed(out_species2, current_species2)
			}

			// update variables at the end of each iteration
			// in the case of ending previous chr, go to start the current chr
			chr_prev_species1 = chr_curr_species1
			chr_prev_species2 = chr_curr_species2
			pos_species1 = 1
			pos_species2 = 1

			current_species1 = bed.Bed{Chrom: chr_curr_species1, ChromStart: pos_species1, ChromEnd: species1_match_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4} // when starting a new chr is when ChromStart needs to be pos_species aka 1, and not the CHromEnd of the previous bed entry
			current_species2 = bed.Bed{Chrom: chr_curr_species2, ChromStart: pos_species2, ChromEnd: species2_match_bed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}

			// before writing bed, test if the bed entries pass filtering criteria
			pass, current_species1.Name, current_species2.Name = gapBedPass(pos_species1, current_species1.ChromStart, current_species1.ChromEnd, pos_species2, current_species2.ChromStart, current_species2.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(out_species1_doNotCalculate, current_species1)
				bed.WriteBed(out_species2_doNotCalculate, current_species2)
			} else {
				bed.WriteBed(out_species1, current_species1)
				bed.WriteBed(out_species2, current_species2)
				// update variables at the end of each iteration
				// only update position if bed passes
				pos_species1 = species1_match_bed[i].ChromEnd
				pos_species2 = species2_match_bed[i].ChromEnd
			}

		} else { // if there is an existing chr

			current_species1 = bed.Bed{Chrom: chr_curr_species1, ChromStart: species1_match_bed[i-1].ChromEnd, ChromEnd: species1_match_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
			current_species2 = bed.Bed{Chrom: chr_curr_species2, ChromStart: species2_match_bed[i-1].ChromEnd, ChromEnd: species2_match_bed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}

			// before writing bed, test if the bed entries pass filtering criteria
			pass, current_species1.Name, current_species2.Name = gapBedPass(pos_species1, current_species1.ChromStart, current_species1.ChromEnd, pos_species2, current_species2.ChromStart, current_species2.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(out_species1_doNotCalculate, current_species1)
				bed.WriteBed(out_species2_doNotCalculate, current_species2)
			} else {
				bed.WriteBed(out_species1, current_species1)
				bed.WriteBed(out_species2, current_species2)
				// update variables at the end of each iteration
				// only update position if bed passes
				pos_species1 = species1_match_bed[i].ChromEnd
				pos_species2 = species2_match_bed[i].ChromEnd
			}
		}
	}

	// after loop, need to write the last entry to output files
	// unless the second to last entry (from the loop) ends at the end of the last chromosome in both species
	// otherwise, write the last entry for both species, so as to keep the number of lines the same between species1 and species2
	if pos_species1 < len(species1_genome_fastaMap[chr_prev_species1]) || pos_species2 < len(species2_genome_fastaMap[chr_prev_species2]) {
		current_species1 = bed.Bed{Chrom: chr_curr_species1, ChromStart: species1_match_bed[len(species1_match_bed)-1].ChromEnd, ChromEnd: len(species1_genome_fastaMap[chr_curr_species1]), Name: "species1_gap", FieldsInitialized: 4}
		current_species2 = bed.Bed{Chrom: chr_curr_species2, ChromStart: species2_match_bed[len(species2_match_bed)-1].ChromEnd, ChromEnd: len(species2_genome_fastaMap[chr_curr_species2]), Name: "species2_gap", FieldsInitialized: 4}
		pass, current_species1.Name, current_species2.Name = gapBedPass(pos_species1, current_species1.ChromStart, current_species1.ChromEnd, pos_species2, current_species2.ChromStart, current_species2.ChromEnd, gapSizeProductLimit)
		if !pass {
			bed.WriteBed(out_species1_doNotCalculate, current_species1)
			bed.WriteBed(out_species2_doNotCalculate, current_species2)
		} else {
			bed.WriteBed(out_species1, current_species1)
			bed.WriteBed(out_species2, current_species2)
		}
	}

	// close output files and check for errors
	err = out_species1.Close()
	exception.PanicOnErr(err)
	err = out_species2.Close()
	exception.PanicOnErr(err)
	err = out_species1_doNotCalculate.Close()
	exception.PanicOnErr(err)
	err = out_species2_doNotCalculate.Close()
	exception.PanicOnErr(err)

	// return output filenames
	return out_species1_filename, out_species2_filename
}

// Step 3: align "gap" sequences
func gapToAlignment(in_species1_gap string, in_species2_gap string, species1_genome string, species2_genome string, species1 string, species2 string, out_filename_prefix string) {
	// read input files
	species1_gap_bed := bed.Read(in_species1_gap)
	species2_gap_bed := bed.Read(in_species2_gap)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_alignment_filename := out_filename_prefix + ".alignment.tsv"
	out_species1_filename := out_filename_prefix + "_" + species1 + "_alignment.bed"
	out_species2_filename := out_filename_prefix + "_" + species2 + "_alignment.bed"
	out_alignment := fileio.EasyCreate(out_alignment_filename)
	out_species1 := fileio.EasyCreate(out_species1_filename)
	out_species2 := fileio.EasyCreate(out_species2_filename)
	var err error

	// initialize variables before loop
	// keep track of sequences
	var species1_seq, species2_seq []dna.Base
	// variables to generate output bed entries
	chr_species1 := ""
	chr_species2 := ""
	pos_species1 := 1
	pos_species2 := 1
	// containers for entries to write to ouput bed files
	var current_species1, current_species2 bed.Bed
	// containers for alignment outputs
	var bestScore int64
	var aln []align.Cigar

	// loop through input gap beds
	for i := range species1_gap_bed {
		// detect for insertions that don't need alignment
		if species1_gap_bed[i].Name == "species1_Insertion" {
			// directly write to output file
			// both species alignment
			bestScore = int64(-600*1 + (-150)*(species1_gap_bed[i].ChromEnd-species1_gap_bed[i].ChromStart-1))                     // without calling alignment function, directly calculate bestScore using gap Open and Extend penalties
			aln = []align.Cigar{{RunLength: int64(species1_gap_bed[i].ChromEnd - species1_gap_bed[i].ChromStart), Op: align.ColD}} // without calling alignment function, directly write alignment cigar
			writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
			// species1 and species2 alignment files
			chr_species1 = species1_gap_bed[i].Chrom
			pos_species1 = species1_gap_bed[i].ChromStart
			bed.WriteBed(out_species1, species1_gap_bed[i])
			pos_species1 += (species1_gap_bed[i].ChromEnd - species1_gap_bed[i].ChromStart)

		} else if species2_gap_bed[i].Name == "species2_Insertion" {
			// directly write to output file
			// both species alignment
			bestScore = int64(-600*1 + (-150)*(species2_gap_bed[i].ChromEnd-species2_gap_bed[i].ChromStart-1))                     // without calling alignment function, directly calculate bestScore using gap Open and Extend penalties
			aln = []align.Cigar{{RunLength: int64(species2_gap_bed[i].ChromEnd - species2_gap_bed[i].ChromStart), Op: align.ColI}} // without calling alignment function, directly write alignment cigar
			writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
			// species1 and species2 alignment files
			chr_species2 = species2_gap_bed[i].Chrom
			pos_species2 = species2_gap_bed[i].ChromStart
			bed.WriteBed(out_species2, species2_gap_bed[i])
			pos_species2 += (species2_gap_bed[i].ChromEnd - species2_gap_bed[i].ChromStart)

		} else {
			// obtain sequences from the genome. To convert bed region (1-based, [closed,open)) to fasta index (0-based, [closed,closed]), subtract 1 from both ChromStart and ChromEnd.
			species1_seq = species1_genome_fastaMap[species1_gap_bed[i].Chrom][(species1_gap_bed[i].ChromStart - 1):(species1_gap_bed[i].ChromEnd - 1)]
			dna.AllToUpper(species1_seq) // convert all bases to uppercase, otherwise get index out of range error in scoring matrix
			species2_seq = species2_genome_fastaMap[species2_gap_bed[i].Chrom][(species2_gap_bed[i].ChromStart - 1):(species2_gap_bed[i].ChromEnd - 1)]
			dna.AllToUpper(species2_seq)

			// align with affineGap, customizeCheckersize, HumanChimpTwoScoreMatrix
			bestScore, aln = align.AffineGap_customizeCheckersize(species1_seq, species2_seq, align.HumanChimpTwoScoreMatrix, -600, -150, 10000, 10000)

			// optional: print results to terminal
			// print alignment score, cigar
			//fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)
			// visualize
			//visualize := align.View(species1_seq, species2_seq, aln)
			//fmt.Println(visualize)

			// write to output file
			// both species alignment
			writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
			// species1 and species2 alignment files
			chr_species1 = species1_gap_bed[i].Chrom
			chr_species2 = species2_gap_bed[i].Chrom
			pos_species1 = species1_gap_bed[i].ChromStart
			pos_species2 = species2_gap_bed[i].ChromStart
			for j := 0; j < len(aln); j++ {
				switch aln[j].Op {
				case align.ColM:
					current_species1 = bed.Bed{Chrom: chr_species1, ChromStart: pos_species1, ChromEnd: pos_species1 + int(aln[j].RunLength), Name: "species1_Match", FieldsInitialized: 4}
					current_species2 = bed.Bed{Chrom: chr_species2, ChromStart: pos_species2, ChromEnd: pos_species2 + int(aln[j].RunLength), Name: "species2_Match", FieldsInitialized: 4}
					bed.WriteBed(out_species1, current_species1)
					bed.WriteBed(out_species2, current_species2)
					pos_species1 += int(aln[j].RunLength)
					pos_species2 += int(aln[j].RunLength)
				case align.ColI:
					current_species2 = bed.Bed{Chrom: chr_species2, ChromStart: pos_species2, ChromEnd: pos_species2 + int(aln[j].RunLength), Name: "species2_Insertion", FieldsInitialized: 4}
					bed.WriteBed(out_species2, current_species2)
					pos_species2 += int(aln[j].RunLength)
				case align.ColD:
					current_species1 = bed.Bed{Chrom: chr_species1, ChromStart: pos_species1, ChromEnd: pos_species1 + int(aln[j].RunLength), Name: "species1_Insertion", FieldsInitialized: 4}
					bed.WriteBed(out_species1, current_species1)
					pos_species1 += int(aln[j].RunLength)
				default:
					log.Fatalf("Unexpected cigar parsing.")
				}
			}
		}
	}

	// close output files and check for errors
	err = out_alignment.Close()
	exception.PanicOnErr(err)
	err = out_species1.Close()
	exception.PanicOnErr(err)
	err = out_species2.Close()
	exception.PanicOnErr(err)
}

// main function: assembles all steps
func globalAlignmentAnchor(in_maf string, species1 string, species2 string, species1_genome string, species2_genome string, gapSizeProductLimit int, chrMap_filename string, out_filename_prefix string, diagonal bool) {
	// process input from out_filename_prefix flag
	// about TrimSuffix: if the suffix string is at the end of the substrate string, then it is trimmed. If the suffix stirng is not at the end of the substrate string, then the substrate string is returned without any change
	if out_filename_prefix == "" {
		out_filename_prefix = strings.TrimSuffix(in_maf, ".maf")
	}

	in_species1_match, in_species2_match := mafToMatch(in_maf, species1, species2, out_filename_prefix, chrMap_filename, diagonal)
	in_species1_gap, in_species2_gap := matchToGap(in_species1_match, in_species2_match, species1_genome, species2_genome, species1, species2, gapSizeProductLimit, out_filename_prefix)
	gapToAlignment(in_species1_gap, in_species2_gap, species1_genome, species2_genome, species1, species2, out_filename_prefix)
}

func usage() {
	fmt.Print(
		"globalAlignmentAnchor - operates on 2 species, takes alignment maf, filters for trusted matches (s lines generated from the same chromosome in both species), and aligns the gap sequences between the trusted matches (affineGap, DefaultScoreMatrix)\n" +
			"in_maf - maf file. Can accept maf describing >2 species, but alignment will be pairwise, aka operating on 2 species\n" +
			"species1, species2 - species names, e.g. hg38. Species1 is target (first line in each maf block); species2 is query (second line in each maf block)\n" +
			"species1_genome, species2_genome - fasta files containing the whole genome of each species. Each fasta sequence is 1 chromosome\n" +
			"chrMap_filename - the name of the file that is the map describing how chromosome names match between species. Each line should be formatted this way: 'species1 chr name'(tab)'species2 chr name'\n" +
			"Usage:\n" +
			"	globalAlignmentAnchor in_maf species1 species2 species1_genome species2_genome chrMap_filename\n" +
			"doNotCalculate flags:\n" +
			"	invalidChromStartOrChromEnd - This bed entry pair is discarded because ChromStart or ChromEnd is invalid\n" +
			"	largeGapSizeMultiple - This bed entry pair is discarded because the gap size of species2 >> the gap size of species1\n" +
			"	large - This bed entry pair is discarded because the product of their sizes is too large\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNum int = 6

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	out_filename_prefix := flag.String("out_filename_prefix", "", "prefix for output filenames")
	diagonal := flag.Bool("diagonal", true, "only allow maf matches close to the chromosome coodrinate diagonal to pass")
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
	chrMap_filename := flag.Arg(5)
	gapSizeProductLimit := 10000000000 // gapSizeProductLimit is currently hardcoded based on align/affineGap tests, 10000000000

	globalAlignmentAnchor(in_maf, species1, species2, species1_genome, species2_genome, gapSizeProductLimit, chrMap_filename, *out_filename_prefix, *diagonal)
}
