package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
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
	defer inputFile.Close()
	for line, doneReading = fileio.EasyNextRealLine(inputFile); !doneReading; line, doneReading = fileio.EasyNextRealLine(inputFile) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
		}
	}
	return seqIdx
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
	//attempt at Dan's len(faOne) suggestion
	//if len(string(faOne)) != 1 || len(string(faTwo)) != 1 {
	//log.Fatalf("multiple sequnces detected in .fa files. This program is designed for .fa files with only 1 sequence in them.\n")
	//}

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
	genomeGraph := simpleGraph.NewGraph()
	//use targetEnd and queryEnd to track position number for each fasta sequence as we move along.
	var targetEnd, queryEnd int = 0, 0
	var curr *simpleGraph.Node
	//TODO: Make sure it can handle non-match at node 0. Fill out annotations for each node (knowing which allele/haplotype).
	//TODO: Handle multi-entry fasta?
	//Creating the first node. This is done independent of all other number because this node has no 'previous' nodes.  All following code leverages the cigar output (second number printed in the Op position of the struct {}) to determine if the alignment returned a 0:match, 1:insertion, or 2:deletion. All indels are relative to the target sequence.
	curr, targetEnd, queryEnd = simpleGraph.FaSeqToNode(faOne, faTwo, targetEnd, queryEnd, aln[0], 0)
	simpleGraph.AddNode(genomeGraph, curr)
	//Drawing the remaining nodes and all edges. Method for adding edges is based on previous nodes.
	for i := 1; i < len(aln); i++ {
		curr, targetEnd, queryEnd = simpleGraph.FaSeqToNode(faOne, faTwo, targetEnd, queryEnd, aln[i], i)
		simpleGraph.AddNode(genomeGraph, curr)
		if aln[i].Op == 0 {
			simpleGraph.AddEdge(genomeGraph.Nodes[i-1], curr, 1)
			simpleGraph.AddEdge(genomeGraph.Nodes[i-2], curr, 0.5)
		} else if aln[i].Op == 1 || aln[i].Op == 2 {
			simpleGraph.AddEdge(genomeGraph.Nodes[i-1], curr, 0.5)
		} else {
			log.Fatalf("Error: cigar.Op = %d is unrecognized...\n", aln[i].Op)
		}
	}
	simpleGraph.PrintGraph(genomeGraph)
}

//raven edited this block to specify only 1 sequnce is expected in each fasta file and add Usage nad options
func usage() {
	fmt.Print(
		"./globalAlignment - chelsea's global alignment\n" +
		" Align 2 .fasta files, each with only 1 sequence\n" +
			"Usage:\n" +
			"	faTarget faQuery\n" +
			"options:\n" +
			"	-faOut=fasta MSA output filename\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNum int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	faOut := flag.String("faOut", "", "name of the MSA output file") //raven added this line
	flag.Parse()

	if len(flag.Args()) != expectedNum {
		flag.Usage()
		log.Fatalf("error, expecting 2 .fasta files to be able to align, only found %d files...\n", len(flag.Args()))
	}

	//read in sequences that should be put in as fasta type files.
	//raven edited this block to save fileio.EasyOpen as file handles, so that the file is only opened 1 time for 2 purposes: faDone and CountSeqIdx
	inputFileOne := fileio.EasyOpen(flag.Arg(0)) //raven's note: EasyOpen returns the type EasyReader
	inputFileTwo := fileio.EasyOpen(flag.Arg(1))
	outFileName := *faOut

	globalAlignment(inputFileOne,inputFileTwo,outFileName)
}
