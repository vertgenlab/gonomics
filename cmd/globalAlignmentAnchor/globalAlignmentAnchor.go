// Command Group: "Linear Alignment Tools"

package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/maf"
)

// helper function: write to tsv file.
func writeToFileHandle(file io.Writer, speciesOne bed.Bed, speciesTwo bed.Bed, score int64, cigar []align.Cigar) {
	var err error
	_, err = fmt.Fprintf(file, "%s\t%s\t%d\t%v\n", speciesOne, speciesTwo, score, cigar)
	exception.PanicOnErr(err)
}

// helper function: make chr name match map.
func makeChrMap(chrMap_filename string) map[string][]string {
	chrMap := make(map[string][]string) // map to hold speciesOne speciesTwo chr name matches
	var chrMapStringSplit []string

	chrMap_string := fileio.Read(chrMap_filename) // read file so each line is 1 string

	for i := range chrMap_string {
		chrMapStringSplit = strings.Split(chrMap_string[i], "\t") // convert each line = 1 string into a slice of 2 strings (assumption: each line only has 2 columns separated by 1 tab). Note that Notepad generates proper tab, but not Atom
		if len(chrMapStringSplit) != 2 {
			log.Fatalf("chrMap did not have 2 columns. Please check that in your chrMap file, each line should be formatted this way: 'speciesOne chr name'(tab)'speciesTwo chr name")
		}

		chrMap[chrMapStringSplit[0]] = append(chrMap[chrMapStringSplit[0]], chrMapStringSplit[1]) // whether or not speciesOne chr name already exists, append to nil or the existing speciesTwo chr name slice the new speciesTwo chr name
	}

	return chrMap
}

// helper function: check if maf entry passes checks.
func matchMafPass(assemblySpeciesOne string, assemblySpeciesTwo string, chromSpeciesOne string, chromSpeciesTwo string, speciesOneSrcSize int, species2_SrcSize int, speciesOneChromStart int, speciesOneChromEnd int, speciesTwoChromStart int, speciesTwoChromEnd int, chrMap map[string][]string, diagonal bool) bool {
	pass := true

	// chrom should match between speciesOne and speciesTwo
	// for each chromSpeciesTwo, go through chrMap to see if it's contained within chromSpeciesOne's matching speciesTwo chr names
	chrMatch := false // separate bool is needed for chrMatch, start out as false
	for _, s := range chrMap[chromSpeciesOne] {
		if chromSpeciesTwo == s {
			chrMatch = true // only become true if chrMatch
		}
	}
	if !chrMatch { // if chrMatch == false, then pass must be set to false
		pass = false
	}

	// maf entry should be roughly diagonal
	// unless user specifies that they do not want the diagonal check
	if diagonal {
		if (float64(speciesTwoChromStart) <= float64(speciesOneChromStart)-0.05*float64(speciesOneSrcSize)) || (float64(speciesTwoChromStart) >= float64(speciesOneChromStart)+0.05*float64(speciesOneSrcSize)) {
			pass = false
		} else if (float64(speciesOneChromStart) <= float64(speciesTwoChromStart)-0.05*float64(species2_SrcSize)) || (float64(speciesOneChromStart) >= float64(speciesTwoChromStart)+0.05*float64(species2_SrcSize)) {
			pass = false
		}
	}

	return pass
}

// helper function: check if gap bed entry passes checks
// need both pos_species (most recently updated alignment position)
// and speciesOneChromStart (the start of the current bed entry).
func gapBedPass(posSpeciesOne int, speciesOneChromStart int, speciesOneChromEnd int, pos_species2 int, speciesTwoChromStart int, speciesTwoChromEnd int, gapSizeProductLimit int) (bool, string, string) {
	pass := true
	speciesOneName := "species1_gap"
	speciesTwoName := "species2_gap"
	// calculate gap size of current bed entry, based on speciesOneChromStart
	speciesOneGapSize := speciesOneChromEnd - speciesOneChromStart
	speciesTwoGapSize := speciesTwoChromEnd - speciesTwoChromStart
	// calculate gap size since the most recently updated alignment position, based on pos_species
	species1_gapSizeBig := speciesOneChromEnd - posSpeciesOne
	species2_gapSizeBig := speciesTwoChromEnd - pos_species2
	// define limits and necessary calculations to evaluate limit
	// gap size product limit should be evaluted based on gapSize
	gapSizeProduct := speciesOneGapSize * speciesTwoGapSize
	// gap size multiple should be evaluated based on gapSizeBig
	gapSizeBigMultiple := 0.0
	if species1_gapSizeBig != 0 {
		gapSizeBigMultiple = float64(species2_gapSizeBig / species1_gapSizeBig)
	}
	gapSizeBigMultipleLimit := 100.00 // gapSizeBigMultipleLimit is currently hardcoded tentatively, 100

	if speciesOneGapSize > 0 && speciesTwoGapSize == 0 { // insertion in speciesOne, no need for alignment
		speciesOneName = "species1_Insertion"
		speciesTwoName = "species2_gap_size0"
	} else if speciesOneGapSize == 0 && speciesTwoGapSize > 0 { // insertion in speciesTwo, no need for alignment
		speciesTwoName = "species2_Insertion"
		speciesOneName = "species1_gap_size0"
	} else if !(speciesOneGapSize > 0 && speciesTwoGapSize > 0) { // need to check, otherwise may have chromEnd<chromStart, leading to gapSize<0. This is one way to make sure in each species esp the non-reference, gap sequence should progress linearly along the chromosome (e.g. alignment match sequence skips around the chromosome, causing gap entries to skip around, ChromStart > ChromEnd)
		pass = false
		speciesOneName = "species1_gap,doNotCalculate_invalidChromStartOrChromEnd"
		speciesTwoName = "species2_gap,doNotCalculate_invalidChromStartOrChromEnd"
	} else if gapSizeBigMultiple > gapSizeBigMultipleLimit { // This is one way to make sure in each species esp the non-reference, gap sequence should progress linearly along the chromosome
		// If the currently aligned sequence in speciesTwo is really far away from the previously aligned sequence in speciesTwo,
		// 100 times further than the currently aligned sequence in speciesOne is from the previously aligned sequence in speciesOne,
		// then the current alignment should be discarded and not "trusted" to be a true alignment,
		// because the currently aligned sequence in speciesTwo is probably a random region further downstream on the chromosome, which breaks synteny.
		// This is only a concern for speciesTwo because speciesOne is always progressing linearly (has synteny) and a big gap in speciesOne might just mean a big insertion in speciesOne.
		pass = false
		speciesOneName = "species1_gap,doNotCalculate_largeGapSizeMultiple"
		speciesTwoName = "species2_gap,doNotCalculate_largeGapSizeMultiple"
		// still check for diagonal. If diagonal, can accept (rescue), but add label
		// this step should come last, after initially accepting some gaps (without aligning), rejecting other gaps, and finally in this step rescuing some rejected gaps to accept them again
		if float64(speciesTwoChromStart) >= 0.95*float64(speciesOneChromStart) && float64(speciesTwoChromStart) <= 1.05*float64(speciesOneChromStart) && float64(speciesTwoChromEnd) >= 0.95*float64(speciesOneChromEnd) && float64(speciesTwoChromEnd) <= 1.05*float64(speciesTwoChromEnd) {
			pass = true
			speciesOneName = "species1_gap_largeGapSize_diagonal"
			speciesTwoName = "species2_gap_largeGapSize_diagonal"
		}
	}

	if gapSizeProduct > gapSizeProductLimit { // the product of the gap sizes needs to be practical for our alignment algorithm. The 2 sequences' product should be <=1E10
		pass = false
		speciesOneName = speciesOneName + ",doNotCalculate_largeGapSizeProduct"
		speciesTwoName = speciesTwoName + ",doNotCalculate_largeGapSizeProduct"
	}

	return pass, speciesOneName, speciesTwoName
}

// Step 1: Filter maf to remove S lines we don't trust, creating filtered maf (aka anchors, or "match")
// not to be confused with cmd/mafFilter, which filters for scores above a threshold.
func mafToMatch(in_maf string, speciesOne string, speciesTwo string, outFilenamePrefix string, chrMap_filename string, diagonal bool) (string, string) {
	mafRecords := maf.Read(in_maf) // read input maf
	chrMap := makeChrMap(chrMap_filename)

	// open output files to write line-by-line and create variable for error
	outMafFilename := outFilenamePrefix + ".filtered.maf"
	out_species1_filename := outFilenamePrefix + "_" + speciesOne + "_match.bed"
	out_species2_filename := outFilenamePrefix + "_" + speciesTwo + "_match.bed"
	outMaf := fileio.EasyCreate(outMafFilename)
	outSpeciesOne := fileio.EasyCreate(out_species1_filename)
	outSpeciesTwo := fileio.EasyCreate(out_species2_filename)
	var err error

	// initialize variables before loop
	// keep track of assembly, chromosome
	var assemblySpeciesOne, assemblySpeciesTwo, chromSpeciesOne, chromSpeciesTwo string
	// containers for entries to write to output files
	var bedSpeciesOne, bedSpeciesTwo bed.Bed

	// loop through input maf
	// I assume pairwise alignment, not >2 species
	// I assume speciesOne is target; speciesTwo is query
	for i := range mafRecords { // each i is a maf block
		// get assembly (e.g. rheMac10), chrom (e.g. chrI)
		assemblySpeciesOne, chromSpeciesOne = maf.SrcToAssemblyAndChrom(mafRecords[i].Species[0].Src)
		// get trusted match coordinates here as well. Output into bed file. Get chrom,start,end,name,score
		bedSpeciesOne = bed.Bed{Chrom: chromSpeciesOne, ChromStart: mafRecords[i].Species[0].SLine.Start, ChromEnd: mafRecords[i].Species[0].SLine.Start + mafRecords[i].Species[0].SLine.Size, Name: "species1_s_filtered_match", Score: int(mafRecords[i].Score), FieldsInitialized: 5}

		for k := 1; k < len(mafRecords[i].Species); k++ { // each k is a line. Start loop at k=1 because that is the lowest possible index to find speciesTwo, which is query
			assemblySpeciesTwo, chromSpeciesTwo = maf.SrcToAssemblyAndChrom(mafRecords[i].Species[k].Src)

			// verify line 0 is indeed speciesOne
			if assemblySpeciesOne != speciesOne {
				log.Fatalf("speciesOne was incorrect. Please check that you have a pairwise maf file, and entered speciesOne and speciesTwo correctly")
			}

			// get s lines
			if mafRecords[i].Species[k].SLine != nil && assemblySpeciesTwo == speciesTwo && mafRecords[i].Species[0].SLine != nil {
				bedSpeciesTwo = bed.Bed{Chrom: chromSpeciesTwo, ChromStart: mafRecords[i].Species[k].SLine.Start, ChromEnd: mafRecords[i].Species[k].SLine.Start + mafRecords[i].Species[k].SLine.Size, Name: "species2_s_filtered_match", Score: int(mafRecords[i].Score), FieldsInitialized: 5}

				// filter out only s lines that we trust to save to filtered maf
				pass := matchMafPass(assemblySpeciesOne, assemblySpeciesTwo, chromSpeciesOne, chromSpeciesTwo, mafRecords[i].Species[0].SLine.SrcSize, mafRecords[i].Species[k].SLine.SrcSize, bedSpeciesOne.ChromStart, bedSpeciesOne.ChromEnd, bedSpeciesTwo.ChromStart, bedSpeciesTwo.ChromEnd, chrMap, diagonal)
				if pass {
					maf.WriteToFileHandle(outMaf, mafRecords[i])
					bed.WriteBed(outSpeciesOne, bedSpeciesOne)
					bed.WriteBed(outSpeciesTwo, bedSpeciesTwo)
				}
			}
		}
	}

	// close output files and check for errors
	err = outMaf.Close()
	exception.PanicOnErr(err)
	err = outSpeciesOne.Close()
	exception.PanicOnErr(err)
	err = outSpeciesTwo.Close()
	exception.PanicOnErr(err)

	// return output filenames
	return out_species1_filename, out_species2_filename
}

// Step 2: Use match to calculate coordinates that still need to be aligned, aka "gap".
func matchToGap(in_species1_match string, in_species2_match string, species1_genome string, species2_genome string, speciesOne string, speciesTwo string, gapSizeProductLimit int, outFilenamePrefix string) (string, string) {
	// read input files
	speciesOneMatchBed := bed.Read(in_species1_match)
	speciesTwoMatchBed := bed.Read(in_species2_match)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_species1_filename := outFilenamePrefix + "_" + speciesOne + "_gap.bed"
	out_species2_filename := outFilenamePrefix + "_" + speciesTwo + "_gap.bed"
	out_species1_doNotCalculate_filename := outFilenamePrefix + "_" + speciesOne + "_gap_doNotCalculate.bed"
	out_species2_doNotCalculate_filename := outFilenamePrefix + "_" + speciesTwo + "_gap_doNotCalculate.bed"
	outSpeciesOne := fileio.EasyCreate(out_species1_filename)
	outSpeciesTwo := fileio.EasyCreate(out_species2_filename)
	outSpeciesOneDoNotCalculate := fileio.EasyCreate(out_species1_doNotCalculate_filename)
	outSpecies2DoNotCalculate := fileio.EasyCreate(out_species2_doNotCalculate_filename)
	var err error

	// initialize variables before loop
	// keep track of chromosome, position
	chr_prev_species1 := speciesOneMatchBed[0].Chrom // initialize prev chr to be the same as curr chr, so the first entry won't have "chrCurrSpeciesOne != chr_prev_species1_species1" and be treated like changing chr
	chrCurrSpeciesOne := speciesOneMatchBed[0].Chrom // initialize curr chr
	chr_prev_species2 := speciesTwoMatchBed[0].Chrom
	chrCurrSpeciesTwo := speciesTwoMatchBed[0].Chrom
	posSpeciesOne := 1 // initialize pos as 1. bed and fa both start at 1
	pos_species2 := 1
	// check for gapBedPass
	var pass bool
	// containers for entries to write to output files
	var currentSpeciesOne, currentSpeciesTwo bed.Bed

	// write i==0 case outside of loop to avoid checking for i==0 in each iteration
	currentSpeciesOne = bed.Bed{Chrom: chrCurrSpeciesOne, ChromStart: posSpeciesOne, ChromEnd: speciesOneMatchBed[0].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
	currentSpeciesTwo = bed.Bed{Chrom: chrCurrSpeciesTwo, ChromStart: pos_species2, ChromEnd: speciesTwoMatchBed[0].ChromStart, Name: "species2_gap", FieldsInitialized: 4}
	// before writing bed, test if the bed entries pass filtering criteria
	pass, currentSpeciesOne.Name, currentSpeciesTwo.Name = gapBedPass(posSpeciesOne, currentSpeciesOne.ChromStart, currentSpeciesOne.ChromEnd, pos_species2, currentSpeciesTwo.ChromStart, currentSpeciesTwo.ChromEnd, gapSizeProductLimit)
	if !pass {
		bed.WriteBed(outSpeciesOneDoNotCalculate, currentSpeciesOne)
		bed.WriteBed(outSpecies2DoNotCalculate, currentSpeciesTwo)
	} else {
		bed.WriteBed(outSpeciesOne, currentSpeciesOne)
		bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
		// update variables at the end of each iteration
		// only update position if bed passes
		posSpeciesOne = speciesOneMatchBed[0].ChromEnd
		pos_species2 = speciesTwoMatchBed[0].ChromEnd
	}

	// starting from i==1
	// loop through input match beds
	// speciesOneMatchBed and speciesTwoMatchBed should have the same number of records, and can be indexed simultaneously
	for i := 1; i < len(speciesOneMatchBed); i++ {
		chrCurrSpeciesOne = speciesOneMatchBed[i].Chrom // set chrCurrSpeciesOne to the new record
		chrCurrSpeciesTwo = speciesTwoMatchBed[i].Chrom

		// calculate the unaligned/gap region before arriving at the aligned/match s line
		if chrCurrSpeciesOne != chr_prev_species1 { // if encounter new chr
			// first finish off the previous chr
			currentSpeciesOne = bed.Bed{Chrom: chr_prev_species1, ChromStart: speciesOneMatchBed[i-1].ChromEnd, ChromEnd: len(species1_genome_fastaMap[chr_prev_species1]), Name: "species1_gap", FieldsInitialized: 4}
			currentSpeciesTwo = bed.Bed{Chrom: chr_prev_species2, ChromStart: speciesTwoMatchBed[i-1].ChromEnd, ChromEnd: len(species2_genome_fastaMap[chr_prev_species2]), Name: "species2_gap", FieldsInitialized: 4}
			pass, currentSpeciesOne.Name, currentSpeciesTwo.Name = gapBedPass(posSpeciesOne, currentSpeciesOne.ChromStart, currentSpeciesOne.ChromEnd, pos_species2, currentSpeciesTwo.ChromStart, currentSpeciesTwo.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(outSpeciesOneDoNotCalculate, currentSpeciesOne)
				bed.WriteBed(outSpecies2DoNotCalculate, currentSpeciesTwo)
			} else {
				bed.WriteBed(outSpeciesOne, currentSpeciesOne)
				bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
			}

			// update variables at the end of each iteration
			// in the case of ending previous chr, go to start the current chr
			chr_prev_species1 = chrCurrSpeciesOne
			chr_prev_species2 = chrCurrSpeciesTwo
			posSpeciesOne = 1
			pos_species2 = 1

			currentSpeciesOne = bed.Bed{Chrom: chrCurrSpeciesOne, ChromStart: posSpeciesOne, ChromEnd: speciesOneMatchBed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4} // when starting a new chr is when ChromStart needs to be pos_species aka 1, and not the CHromEnd of the previous bed entry
			currentSpeciesTwo = bed.Bed{Chrom: chrCurrSpeciesTwo, ChromStart: pos_species2, ChromEnd: speciesTwoMatchBed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}

			// before writing bed, test if the bed entries pass filtering criteria
			pass, currentSpeciesOne.Name, currentSpeciesTwo.Name = gapBedPass(posSpeciesOne, currentSpeciesOne.ChromStart, currentSpeciesOne.ChromEnd, pos_species2, currentSpeciesTwo.ChromStart, currentSpeciesTwo.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(outSpeciesOneDoNotCalculate, currentSpeciesOne)
				bed.WriteBed(outSpecies2DoNotCalculate, currentSpeciesTwo)
			} else {
				bed.WriteBed(outSpeciesOne, currentSpeciesOne)
				bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
				// update variables at the end of each iteration
				// only update position if bed passes
				posSpeciesOne = speciesOneMatchBed[i].ChromEnd
				pos_species2 = speciesTwoMatchBed[i].ChromEnd
			}
		} else { // if there is an existing chr
			currentSpeciesOne = bed.Bed{Chrom: chrCurrSpeciesOne, ChromStart: speciesOneMatchBed[i-1].ChromEnd, ChromEnd: speciesOneMatchBed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
			currentSpeciesTwo = bed.Bed{Chrom: chrCurrSpeciesTwo, ChromStart: speciesTwoMatchBed[i-1].ChromEnd, ChromEnd: speciesTwoMatchBed[i].ChromStart, Name: "species2_gap", FieldsInitialized: 4}

			// before writing bed, test if the bed entries pass filtering criteria
			pass, currentSpeciesOne.Name, currentSpeciesTwo.Name = gapBedPass(posSpeciesOne, currentSpeciesOne.ChromStart, currentSpeciesOne.ChromEnd, pos_species2, currentSpeciesTwo.ChromStart, currentSpeciesTwo.ChromEnd, gapSizeProductLimit)
			if !pass {
				bed.WriteBed(outSpeciesOneDoNotCalculate, currentSpeciesOne)
				bed.WriteBed(outSpecies2DoNotCalculate, currentSpeciesTwo)
			} else {
				bed.WriteBed(outSpeciesOne, currentSpeciesOne)
				bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
				// update variables at the end of each iteration
				// only update position if bed passes
				posSpeciesOne = speciesOneMatchBed[i].ChromEnd
				pos_species2 = speciesTwoMatchBed[i].ChromEnd
			}
		}
	}

	// after loop, need to write the last entry to output files
	// unless the second to last entry (from the loop) ends at the end of the last chromosome in both species
	// otherwise, write the last entry for both species, so as to keep the number of lines the same between speciesOne and speciesTwo
	if posSpeciesOne < len(species1_genome_fastaMap[chr_prev_species1]) || pos_species2 < len(species2_genome_fastaMap[chr_prev_species2]) {
		currentSpeciesOne = bed.Bed{Chrom: chrCurrSpeciesOne, ChromStart: speciesOneMatchBed[len(speciesOneMatchBed)-1].ChromEnd, ChromEnd: len(species1_genome_fastaMap[chrCurrSpeciesOne]), Name: "species1_gap", FieldsInitialized: 4}
		currentSpeciesTwo = bed.Bed{Chrom: chrCurrSpeciesTwo, ChromStart: speciesTwoMatchBed[len(speciesTwoMatchBed)-1].ChromEnd, ChromEnd: len(species2_genome_fastaMap[chrCurrSpeciesTwo]), Name: "species2_gap", FieldsInitialized: 4}
		pass, currentSpeciesOne.Name, currentSpeciesTwo.Name = gapBedPass(posSpeciesOne, currentSpeciesOne.ChromStart, currentSpeciesOne.ChromEnd, pos_species2, currentSpeciesTwo.ChromStart, currentSpeciesTwo.ChromEnd, gapSizeProductLimit)
		if !pass {
			bed.WriteBed(outSpeciesOneDoNotCalculate, currentSpeciesOne)
			bed.WriteBed(outSpecies2DoNotCalculate, currentSpeciesTwo)
		} else {
			bed.WriteBed(outSpeciesOne, currentSpeciesOne)
			bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
		}
	}

	// close output files and check for errors
	err = outSpeciesOne.Close()
	exception.PanicOnErr(err)
	err = outSpeciesTwo.Close()
	exception.PanicOnErr(err)
	err = outSpeciesOneDoNotCalculate.Close()
	exception.PanicOnErr(err)
	err = outSpecies2DoNotCalculate.Close()
	exception.PanicOnErr(err)

	// return output filenames
	return out_species1_filename, out_species2_filename
}

// Step 3: align "gap" sequences.
func gapToAlignment(in_species1_gap string, in_species2_gap string, species1_genome string, species2_genome string, speciesOne string, speciesTwo string, outFilenamePrefix string) {
	// read input files
	species1_gap_bed := bed.Read(in_species1_gap)
	species2_gap_bed := bed.Read(in_species2_gap)
	species1_genome_fa := fasta.Read(species1_genome)
	species1_genome_fastaMap := fasta.ToMap(species1_genome_fa)
	species2_genome_fa := fasta.Read(species2_genome)
	species2_genome_fastaMap := fasta.ToMap(species2_genome_fa)

	// open output files to write line-by-line and create variable for error
	out_alignment_filename := outFilenamePrefix + ".alignment.tsv"
	out_species1_filename := outFilenamePrefix + "_" + speciesOne + "_alignment.bed"
	out_species2_filename := outFilenamePrefix + "_" + speciesTwo + "_alignment.bed"
	out_alignment := fileio.EasyCreate(out_alignment_filename)
	outSpeciesOne := fileio.EasyCreate(out_species1_filename)
	outSpeciesTwo := fileio.EasyCreate(out_species2_filename)
	var err error

	// initialize variables before loop
	// keep track of sequences
	var species1_seq, species2_seq []dna.Base
	// variables to generate output bed entries
	chr_species1 := ""
	chr_species2 := ""
	posSpeciesOne := 1
	pos_species2 := 1
	// containers for entries to write to output bed files
	var currentSpeciesOne, currentSpeciesTwo bed.Bed
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
			// speciesOne and speciesTwo alignment files
			chr_species1 = species1_gap_bed[i].Chrom
			posSpeciesOne = species1_gap_bed[i].ChromStart
			bed.WriteBed(outSpeciesOne, species1_gap_bed[i])
			posSpeciesOne += (species1_gap_bed[i].ChromEnd - species1_gap_bed[i].ChromStart)
		} else if species2_gap_bed[i].Name == "species2_Insertion" {
			// directly write to output file
			// both species alignment
			bestScore = int64(-600*1 + (-150)*(species2_gap_bed[i].ChromEnd-species2_gap_bed[i].ChromStart-1))                     // without calling alignment function, directly calculate bestScore using gap Open and Extend penalties
			aln = []align.Cigar{{RunLength: int64(species2_gap_bed[i].ChromEnd - species2_gap_bed[i].ChromStart), Op: align.ColI}} // without calling alignment function, directly write alignment cigar
			writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
			// speciesOne and speciesTwo alignment files
			chr_species2 = species2_gap_bed[i].Chrom
			pos_species2 = species2_gap_bed[i].ChromStart
			bed.WriteBed(outSpeciesTwo, species2_gap_bed[i])
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
			// fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)
			// visualize
			// visualize := align.View(species1_seq, species2_seq, aln)
			// fmt.Println(visualize)

			// write to output file
			// both species alignment
			writeToFileHandle(out_alignment, species1_gap_bed[i], species2_gap_bed[i], bestScore, aln)
			// speciesOne and speciesTwo alignment files
			chr_species1 = species1_gap_bed[i].Chrom
			chr_species2 = species2_gap_bed[i].Chrom
			posSpeciesOne = species1_gap_bed[i].ChromStart
			pos_species2 = species2_gap_bed[i].ChromStart
			for j := 0; j < len(aln); j++ {
				switch aln[j].Op {
				case align.ColM:
					currentSpeciesOne = bed.Bed{Chrom: chr_species1, ChromStart: posSpeciesOne, ChromEnd: posSpeciesOne + int(aln[j].RunLength), Name: "species1_Match", FieldsInitialized: 4}
					currentSpeciesTwo = bed.Bed{Chrom: chr_species2, ChromStart: pos_species2, ChromEnd: pos_species2 + int(aln[j].RunLength), Name: "species2_Match", FieldsInitialized: 4}
					bed.WriteBed(outSpeciesOne, currentSpeciesOne)
					bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
					posSpeciesOne += int(aln[j].RunLength)
					pos_species2 += int(aln[j].RunLength)
				case align.ColI:
					currentSpeciesTwo = bed.Bed{Chrom: chr_species2, ChromStart: pos_species2, ChromEnd: pos_species2 + int(aln[j].RunLength), Name: "species2_Insertion", FieldsInitialized: 4}
					bed.WriteBed(outSpeciesTwo, currentSpeciesTwo)
					pos_species2 += int(aln[j].RunLength)
				case align.ColD:
					currentSpeciesOne = bed.Bed{Chrom: chr_species1, ChromStart: posSpeciesOne, ChromEnd: posSpeciesOne + int(aln[j].RunLength), Name: "species1_Insertion", FieldsInitialized: 4}
					bed.WriteBed(outSpeciesOne, currentSpeciesOne)
					posSpeciesOne += int(aln[j].RunLength)
				default:
					log.Fatalf("Unexpected cigar parsing.")
				}
			}
		}
	}

	// close output files and check for errors
	err = out_alignment.Close()
	exception.PanicOnErr(err)
	err = outSpeciesOne.Close()
	exception.PanicOnErr(err)
	err = outSpeciesTwo.Close()
	exception.PanicOnErr(err)
}

// main function: assembles all steps.
func globalAlignmentAnchor(in_maf string, speciesOne string, speciesTwo string, species1_genome string, species2_genome string, gapSizeProductLimit int, chrMap_filename string, outFilenamePrefix string, diagonal bool) {
	// process input from outFilenamePrefix flag
	// about TrimSuffix: if the suffix string is at the end of the substrate string, then it is trimmed. If the suffix stirng is not at the end of the substrate string, then the substrate string is returned without any change
	if outFilenamePrefix == "" {
		outFilenamePrefix = strings.TrimSuffix(in_maf, ".maf")
	}

	in_species1_match, in_species2_match := mafToMatch(in_maf, speciesOne, speciesTwo, outFilenamePrefix, chrMap_filename, diagonal)
	in_species1_gap, in_species2_gap := matchToGap(in_species1_match, in_species2_match, species1_genome, species2_genome, speciesOne, speciesTwo, gapSizeProductLimit, outFilenamePrefix)
	gapToAlignment(in_species1_gap, in_species2_gap, species1_genome, species2_genome, speciesOne, speciesTwo, outFilenamePrefix)
}

func usage() {
	fmt.Print(
		"globalAlignmentAnchor - operates on 2 species, takes alignment maf, filters for trusted matches (s lines generated from the same chromosome in both species), and aligns the gap sequences between the trusted matches (affineGap, DefaultScoreMatrix)\n" +
			"in_maf - maf file. Can accept maf describing >2 species, but alignment will be pairwise, aka operating on 2 species\n" +
			"speciesOne, speciesTwo - species names, e.g. hg38. Species1 is target (first line in each maf block); speciesTwo is query (second line in each maf block)\n" +
			"species1_genome, species2_genome - fasta files containing the whole genome of each species. Each fasta sequence is 1 chromosome\n" +
			"chrMap_filename - the name of the file that is the map describing how chromosome names match between species. Each line should be formatted this way: 'speciesOne chr name'(tab)'speciesTwo chr name'\n" +
			"Usage:\n" +
			"	globalAlignmentAnchor in_maf speciesOne speciesTwo species1_genome species2_genome chrMap_filename\n" +
			"doNotCalculate flags:\n" +
			"	invalidChromStartOrChromEnd - This bed entry pair is discarded because ChromStart or ChromEnd is invalid\n" +
			"	largeGapSizeMultiple - This bed entry pair is discarded because the gap size of speciesTwo >> the gap size of speciesOne\n" +
			"	large - This bed entry pair is discarded because the product of their sizes is too large\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNum int = 6

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	outFilenamePrefix := flag.String("outFilenamePrefix", "", "prefix for output filenames")
	diagonal := flag.Bool("diagonal", true, "only allow maf matches close to the chromosome coodrinate diagonal to pass")
	flag.Parse()

	if len(flag.Args()) != expectedNum {
		flag.Usage()
		log.Fatalf("error, expecting 5 command line arguments, only found %d\n", len(flag.Args()))
	}

	in_maf := flag.Arg(0)
	speciesOne := flag.Arg(1)
	speciesTwo := flag.Arg(2)
	species1_genome := flag.Arg(3)
	species2_genome := flag.Arg(4)
	chrMap_filename := flag.Arg(5)
	gapSizeProductLimit := 10000000000 // gapSizeProductLimit is currently hardcoded based on align/affineGap tests, 10000000000

	globalAlignmentAnchor(in_maf, speciesOne, speciesTwo, species1_genome, species2_genome, gapSizeProductLimit, chrMap_filename, *outFilenamePrefix, *diagonal)
}
