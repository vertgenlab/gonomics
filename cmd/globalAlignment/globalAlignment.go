package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"strings" //raven added this line for CountSeqIdx
	"os" //raven added this line for MSA Fasta output
)

//raven edited this block to specify only 1 sequnce is expected in each fasta file
func usage() {
	fmt.Print(
	"./globalAlignment - chelsea's global alignment\n" +
	"Usage:\n" +
	" Align 2 .fasta files, each with only 1 sequence\n" +
	"options:\n")
	flag.PrintDefaults()
}

//raven wrote this block to count sequences based on the Read function in gonomics/fasta/fasta.go
//maybe there are existing functions that do this already? I didn't find an existing function that returns seqIdx
func CountSeqIdx(filename string) int {
	var line string
	var seqIdx int = 0 //I know in Read seqIdx is int64 and starts with -1, but I am using it differently here
	var doneReading bool = false
	file := fileio.EasyOpen(filename)
	defer file.Close()
	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		if strings.HasPrefix(line, ">") {
			seqIdx++
		}
	}
	return seqIdx
}

func main() {
	faOut := flag.String("faOut","globalAlignment_result.fa","name of the MSA output file") //raven added this line
	flag.Parse()
	var expectedNum int = 2
	if len(flag.Args()) != expectedNum {
		log.Fatalf("error, expecting 2 .fasta files to be able to align, only found %d files...\n", len(flag.Args()))
	}

	//read in sequences that should be put in as fasta type files.
	faOne, faDoneOne := fasta.NextFasta(fileio.EasyOpen(flag.Arg(0)))
	faTwo, faDoneTwo := fasta.NextFasta(fileio.EasyOpen(flag.Arg(1)))

	//fmt.Printf("%v \n %v \n", faOne, faTwo)
	if faDoneOne || faDoneTwo {
		log.Fatalf("error, unable to read .fa files, check for > symbol before each name and that each fasta entry is only one line. Check to make sure there are only four bases or N, globalAlignment is unable to use anything outside A-T-G-C-N.\n")
	}

	//raven added this block
	numSeqOne := CountSeqIdx(flag.Arg(0))
	numSeqTwo := CountSeqIdx(flag.Arg(1))
	if numSeqOne > 1 || numSeqTwo > 1 {
		log.Fatalf("multiple sequnces detected in .fa files: %v sequences in the first .fa file and %v sequences in the second .fa file. This program is designed for .fa files with only 1 sequence in them\n",numSeqOne,numSeqTwo)
	}

	//needleman wunsch (global alignment)
	bestScore, aln := align.ConstGap(faOne.Seq, faTwo.Seq, align.HumanChimpTwoScoreMatrix, -430)
	fmt.Printf("Alignment score is %d, cigar is %v \n", bestScore, aln)

	//visualize
	visualize := align.View(faOne.Seq, faTwo.Seq, aln)
	fmt.Println(visualize)

	//raven added this block to put visualized alignment output into MSA Fasta
	visualizeOutput := ">" + faOne.Name + "\n" + strings.Split(visualize,"\n")[0] + "\n" + ">" + faTwo.Name + "\n" + strings.Split(visualize,"\n")[1]
	outFileName := *faOut
	outFile,err1 := os.Create(outFileName) //fileio.EasyCreate and other fileio tools should work, but I don't really know how to use them
	if err1 != nil {
		log.Fatalf("Write to file failed on step 1\n")
	}
	defer outFile.Close()
	_,err2 := outFile.WriteString(visualizeOutput)
	if err2 != nil {
		log.Fatalf("Write to file failed on step 2\n")
	}

	//graph
	genomeGraph := cigarToGraph(faOne, faTwo, aln)
	simpleGraph.PrintGraph(genomeGraph)

}

//faOne target, faTwo query
//cigar to graph takes two fasta files as an input, performs a global alignment, and produces genome graph
func cigarToGraph(target *fasta.Fasta, query *fasta.Fasta, aln []align.Cigar) *simpleGraph.SimpleGraph {
	answer := simpleGraph.NewGraph()
	//use targetEnd and queryEnd to track position number for each fasta sequence as we move along.
	var targetEnd, queryEnd int = 0, 0
	var curr *simpleGraph.Node
	//TODO: Make sure it can handle non-match at node 0. Fill out annotations for each node (knowing which allele/haplotype).
	//TODO: Handle multi-entry fasta?

	//Creating the first node. This is done independent of all other number because this node has no 'previous' nodes.  All following code leverages the cigar output (second number printed in the Op position of the struct {}) to determine if the alignment returned a 0:match, 1:insertion, or 2:deletion. All indels are relative to the target sequence.
	curr, targetEnd, queryEnd = simpleGraph.FaSeqToNode(target, query, targetEnd, queryEnd, aln[0], 0)
	simpleGraph.AddNode(answer, curr)
	//Drawing the remaining nodes and all edges. Method for adding edges is based on previous nodes.
	for i := 1; i < len(aln); i++ {
		curr, targetEnd, queryEnd = simpleGraph.FaSeqToNode(target, query, targetEnd, queryEnd, aln[i], i)
		simpleGraph.AddNode(answer, curr)
		if aln[i].Op == 0 {
			simpleGraph.AddEdge(answer.Nodes[i-1], curr, 1)
			simpleGraph.AddEdge(answer.Nodes[i-2], curr, 0.5)
		} else if aln[i].Op == 1 || aln[i].Op == 2 {
			simpleGraph.AddEdge(answer.Nodes[i-1], curr, 0.5)
		} else {
			log.Fatalf("Error: cigar.Op = %d is unrecognized...\n", aln[i].Op)
		}
	}
	return answer
}
