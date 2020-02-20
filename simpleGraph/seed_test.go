package simpleGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"testing"
)

//TODO: consider dot format to graph
//     --> node_1            --> node_4
//node_0          --> node_3             --> node_6
//     --> node_2           -->  node_5

func TestSeedExtend(t *testing.T) {
	gg := NewGraph()
	newNode := &Node{Id: 0, Seq: dna.StringToBases("CCCATTCTTCCCTCCACAAGAT")}
	AddNode(gg, newNode)
	snp1 := &Node{Id: 1, Seq: dna.StringToBases("C")}
	AddNode(gg, snp1)
	AddEdge(newNode, snp1, 0.5)
	snp2 := &Node{Id: 2, Seq: dna.StringToBases("C")}
	AddNode(gg, snp2)
	AddEdge(newNode, snp2, 0.5)
	match2 := &Node{Id: 3, Seq: dna.StringToBases("CAGATCCCCGTGTTTGGACCCCATAAGAATGAGAAGA")}
	AddNode(gg, match2)
	AddEdge(snp1, match2, 0.5)
	AddEdge(snp2, match2, 0.5)
	snp3 := &Node{Id: 4, Seq: dna.StringToBases("AT")}
	AddNode(gg, snp3)
	AddEdge(match2, snp3, 0.5)
	snp4 := &Node{Id: 5, Seq: dna.StringToBases("AT")}
	AddNode(gg, snp4)
	AddEdge(match2, snp4, 0.5)
	lastMatch := &Node{Id: 6, Seq: dna.StringToBases("CCTCTCCGGTTTTGTTATCTGTTCTGAATGTGAAACTGGGATATCCACGTTGATCTCGGCTAAAATATAACAAACCTAAACAATTTGA")}
	AddNode(gg, lastMatch)
	AddEdge(snp3, lastMatch, 0.5)
	AddEdge(snp4, lastMatch, 0.5)
	var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("CCCATTCTTCCCTCCACAAGATCCAGATCCCCGTGTTTGGACCCCATAAGAATGAGAAGAATCCTCTCCGGTTTTGTTATCTGTTCTGAATGTGAAACTGGGATATCCACGTTGATCTCGGCTAAAATATAACAAACCTAAACAATTTGA"), Qual: []rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ")}
	seed := &SeedDev{TargetId: 3, TargetStart: 0, QueryStart: 23, Length: 10}
	hits := extendSeedTogether(seed, gg, hippo)
	printSeedDev(hits)
	log.Printf("Number of expanded seeds is: %d", len(hits))
	seq1 := dna.StringToBases("TCCACCTAGAGGAGCCTGTTCTAGAACCGATAACCCCCGTTCAACCTCACCTCCCCTTGTTAATACCGCCTATATACCACCGTCGTCAGCTTACCCTGTGAGGGACTAATAGTAAGCTAAACTGGTATAACCCTAAACGTCAGGTCGAGGTGTAGCGTATGTGGAGGGAAGAAATGGGCTACATTC-GCTACAAATAGCGAACACGAATGATGTCCTGAAATGTACATCTGAAGGAGGATTTAGCAGTAAGTAGAAAATAGAGTGT")
	seq2 := dna.StringToBases("TCCACCTAGAGGAGCCTGTTCTAGAACCGATAACCACCGTTCAACCTTACCTCCCCTTGTTAATACCGCCTATATACCACCGTCGTCAGCTTACCCTGTGAGGGACTAATAGTAAGCTAAACT-------------------GGTCGAGGTGTAGCGAATGAGGAGGGAAGAAATGGGCTACATTTGGCTAC-AGCAGCGAACACGAATGATGTCTTGAAACGTACATCTGAAGGAGGATTTAGCAGTAAGTACAAAATAGAGTGT")
	log.Printf("Length of 2 sequences are: %d, %d\n", len(seq1), len(seq2))
}
