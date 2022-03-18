package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception" //raven added this line for file Close, exception.PanicOnErr
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"log"
	"os"      //raven added this line for MSA Fasta output
	"strings" //raven added this line for CountSeqIdx
)

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
func GlobalAlignment_CigarToBed(inputFileOne *fileio.EasyReader, inputFileTwo *fileio.EasyReader, outFa string, outIns_bed string, outDel_bed string, FirstPos_InsBed int, FirstPos_DelBed int, Chrom string) {
	//make sure files meet the usage requirements of globalAlignment.go
	faOne, faDoneOne := fasta.NextFasta(inputFileOne)
	faTwo, faDoneTwo := fasta.NextFasta(inputFileTwo)
	fasta.ToUpper(faOne) //raven added these 2 lines to prevent lowercase letters from indexing out of range in align/linearGap.go
	fasta.ToUpper(faTwo)
	if faDoneOne || faDoneTwo {
		log.Fatalf("error, unable to read .fa files, check for > symbol before each name and that each fasta entry is only one line. Check to make sure there are only four bases or N, globalAlignment is unable to use anything outside A-T-G-C-N.\n")
	}
	numSeqOne := CountSeqIdx(inputFileOne)
	numSeqTwo := CountSeqIdx(inputFileTwo)
	if numSeqOne > 1 || numSeqTwo > 1 {
		log.Fatalf("multiple sequnces detected in .fa files: %v sequences in the first .fa file and %v sequences in the second .fa file. This program is designed for .fa files with only 1 sequence in them\n", numSeqOne, numSeqTwo)
	}

	//needleman wunsch (global alignment)
	//raven's note: scoring matrices are in align/align.go, e.g. HumanChimp matrix, same as what UCSC Genome Browser uses
	//raven's note: ConstGap is a function in align/linearGap.go, with only constant gapPen (gap penalty)
	//bestScore, aln := align.ConstGap(faOne.Seq, faTwo.Seq, align.HumanChimpTwoScoreMatrix, -430)
	//fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)
	//raven's edit: try AffineGap in align/affineGap.go, with diffrent gapOpen and gapExtend penalty, since gapExtend<gapOpen, encourages big gap in order to align better in unbroken chunks
	bestScore, aln := align.AffineGap(faOne.Seq, faTwo.Seq, align.HumanChimpTwoScoreMatrix, -600, -150)
	fmt.Printf("Using AffineGap, Alignment score is %d, cigar is %v \n", bestScore, aln)

	//raven's cigarToBed main body starts here
	//insertion bed
	insBed := fileio.EasyCreate(outIns_bed)
	//TODO: now have to write FirstPos, but find a way to grab from fasta header
	ChromCurrent := FirstPos_InsBed - 1 //initialize a variable to keep track of current position on chromosome
	ChromStart := FirstPos_InsBed - 1   //initialize variables for ChromStart and ChromEnd, will need to set them equal to one another during updates, but may be more efficient to have fewer variables
	ChromEnd := FirstPos_InsBed - 1
	BedEntry := bed.Bed{Chrom: Chrom, ChromStart: ChromStart, ChromEnd: ChromEnd, Name: "ins", FieldsInitialized: 4} //TODO: now just write "chr1", but find a way to grab from fasta header?
	for i := 0; i < len(aln)-1; i++ {                                                                                //loop through each entry in the []cigar, start from 0, end at len-2 so can still assess aln[i+1]
		if aln[i].Op == align.ColM && aln[i+1].Op == align.ColI { //if there is an M before I, then write to bed
			ChromStart = ChromCurrent + int(aln[i].RunLength) + 1 //RunLengths are int64, but bed fields like ChromStart are int, so need to convert RunLengths to int, get the position at which I starts
			ChromEnd = ChromStart + int(aln[i+1].RunLength)       //get the last position that is I
			BedEntry = bed.Bed{Chrom: Chrom, ChromStart: ChromStart, ChromEnd: ChromEnd, Name: "ins", FieldsInitialized: 4}
			bed.WriteBed(insBed, BedEntry)
		}
		if aln[i].Op != align.ColD { //in insertion bed, only need to update ChromCurrent if the cigar fragment is M or I, not D
			ChromCurrent += int(aln[i].RunLength)
		}
	}
	err := insBed.Close()
	exception.PanicOnErr(err)

	//deletion bed
	delBed := fileio.EasyCreate(outDel_bed) //TODO: add flag for output file name
	//TODO: now have to write FirstPos, but find a way to grab from fasta header
	ChromCurrent = FirstPos_DelBed - 1 //reset the same 3 variables for position
	ChromStart = FirstPos_DelBed - 1
	ChromEnd = FirstPos_DelBed - 1
	for i := 0; i < len(aln)-1; i++ {
		if aln[i].Op == align.ColM && aln[i+1].Op == align.ColI {
			ChromStart = ChromCurrent + int(aln[i].RunLength) //the position before I starts
			ChromEnd = ChromStart + 1                         //add 1 so that the length of the bed entry is 1
			BedEntry = bed.Bed{Chrom: Chrom, ChromStart: ChromStart, ChromEnd: ChromEnd, Name: "del", FieldsInitialized: 4}
			bed.WriteBed(delBed, BedEntry)
		}
		if aln[i].Op != align.ColI { //in deletion bed, only update ChromCurrent if the cigar fragment is M or D, not I
			ChromCurrent += int(aln[i].RunLength)
		}
	}
	err = delBed.Close()
	exception.PanicOnErr(err)

	//visualize
	visualize := align.View(faOne.Seq, faTwo.Seq, aln)
	fmt.Println(visualize)

	//raven added this block to put visualized alignment output into MSA Fasta
	if outFa != "" { //only write to file if a filename is given in the faOut command line option
		outFaFile, err := os.Create(outFa) //fileio.EasyCreate and other fileio tools should work, but I don't really know how to use them
		if err != nil {
			log.Fatalf("Write to file failed on step 1\n")
		}
		visualizeOutput := ">" + faOne.Name + "\n" + strings.Split(visualize, "\n")[0] + "\n" + ">" + faTwo.Name + "\n" + strings.Split(visualize, "\n")[1] + "\n"
		_, err = outFaFile.WriteString(visualizeOutput)
		if err != nil {
			log.Fatalf("Write to file failed on step 2\n")
		}
		outFaFile.Close() //commented out defer close
	}

	//cigar to graph
	//genomeGraph := cigarToGraph(faOne, faTwo, aln)
	//genomeGraph.PrintGraph(genomeGraph)
}

//raven edited this block to specify only 1 sequnce is expected in each fasta file and add Usage nad options
func usage() {
	fmt.Print(
		"./cigarToBed\n" +
			" Uses globalAlignment, affineGap (instead of constGap) to align 2 .fasta files, each with only 1 sequence, then convert cigars to ins and del beds\n" +
			"Usage:\n" +
			" cigartobed target.fasta query.fasta\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNum int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	faOut := flag.String("faOut", "", "fasta MSA output filename")
	insBedOut := flag.String("insBedOut", "ins.bed", "insertion bed output filename")
	delBedOut := flag.String("delBedOut", "del.bed", "deletion bed output filename")
	FirstPos_Ins := flag.Int("FirstPos_Ins", 1, "first base position of the query sequence, which will be used to make the insertion bed")
	FirstPos_Del := flag.Int("FirstPos_Del", 1, "first base position of the target sequence, which will be used to make the deletion bed")
	Chr := flag.String("Chr", "chr1", "chromosome name")
	flag.Parse()

	if len(flag.Args()) != expectedNum {
		flag.Usage()
		log.Fatalf("error, expecting 2 .fasta files to be able to align, only found %d files...\n", len(flag.Args()))
	}

	//read in sequences that should be put in as fasta type files.
	//raven edited this block to save fileio.EasyOpen as file handles, so that the file is only opened 1 time for 2 purposes: faDone and CountSeqIdx
	inputFileOne := fileio.EasyOpen(flag.Arg(0)) //raven's note: EasyOpen returns the type EasyReader
	inputFileTwo := fileio.EasyOpen(flag.Arg(1))

	GlobalAlignment_CigarToBed(inputFileOne, inputFileTwo, *faOut, *insBedOut, *delBedOut, *FirstPos_Ins, *FirstPos_Del, *Chr)
}
