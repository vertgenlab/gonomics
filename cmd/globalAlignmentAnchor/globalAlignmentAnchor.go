// Command Group: "Linear Alignment Tools"

package main

import (
	"flag"
	"fmt"
	//"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/maf"
	"github.com/vertgenlab/gonomics/bed"
	"log"
	//"os"
	"strings"
)

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
	species1_bed := bed.Read(species1_match_bed_filename)
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
	for i := range species1_bed {
		chr_curr = species1_bed[i].Chrom //set chr_curr to the new record
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
		current_species1 := bed.Bed{Chrom: chr_curr, ChromStart: pos_species1, ChromEnd: species1_bed[i].ChromStart, Name: "species1_gap", FieldsInitialized: 4}
		bed.WriteBed(out_species1.File, current_species1)
		bed.WriteBed(out_species2.File, current_species2)

		//update variables at the end of each iteration
		pos_species1 = species1_bed[i].ChromEnd
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

/*
//Step 3: globalAlignment lowMem for non-anchor sequences
//raven did not put this helper function into the globalAlignment function because it is used twice within the globalAlignment function
//raven wrote this block to count sequences based on the Read function in gonomics/fasta/fasta.go
//raven changed the input variable from filename string to inputFile EasyReader, so that the file is only opened 1 time for 2 purposes: faDone and CountSeqIdx
//resources
	species1_genome_fastaMap[i].Name
	records[i].Seq[start:end] // I believe Seq starts index at 0, according to fasta.go WriteFasta function
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

	//Creating the first node. This is done independent of all other number because this node has no 'previous' nodes.  All following code leverages the cigar output (second number printed in the Op position of the struct {}) to determine if the alignment returned a 0:match, 1:species1ertion, or 2:species2etion. All inspecies2s are relative to the target sequence.
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
*/

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
}
