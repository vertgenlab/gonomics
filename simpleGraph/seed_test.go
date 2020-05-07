package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"testing"
)

//TODO: consider dot format to graph
//     --> node_1            --> node_4
//node_0          --> node_3             --> node_6
//     --> node_2           -->  node_5
func TestScoreSeed(t *testing.T) {
	fq := fastq.Fastq{Name: "Test", Seq: dna.StringToBases("AAAAGCTCAGCTCCCATCGAGCTGAAAATGTTCCCCCCCACAGTTTTATCATGCTGTAAATGAAGGTGGGATGTCTCAACAAATAGTTCCCAGCCTTCTTTTTCTCTTTTGTGAAAAAGGCCAAAGGAATTCAGACTGTCACGTAATTTA")}
	seedScore := perfectMatch(&fq, HumanChimpTwoScoreMatrix)
	log.Printf("Score of this read is: %d", seedScore)
}
func TestSeed(t *testing.T) {
	gg := NewGraph()
	newNode := &Node{Id: 0, Seq: dna.StringToBases("CCCATTCTTCCCTCCACAAGAT")}
	AddNode(gg, newNode)
	snp1 := &Node{Id: 1, Seq: dna.StringToBases("C")}
	AddNode(gg, snp1)
	AddEdge(newNode, snp1, 0.5)
	snp2 := &Node{Id: 2, Seq: dna.StringToBases("T")}
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
	var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("CCCATTCTTCCCTCCACAAGATCCAGATCCCCGTGTTTGGACCCCATAAGAATGAGAAGAATCCTCTCCGGTTTTGTTATCTGTTCTGAATGTGAAACTGGGATATCCACGTTGATCTCGGCTAAAATATAACAAACCTAAACAATTTGA"), Qual: fastq.ToQualUint8([]rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"))}
	seed := &SeedDev{TargetId: 3, TargetStart: 0, QueryStart: 23, Length: 10}
	hits := extendSeedTogether(seed, gg, hippo)
	printSeedDev(hits)

	log.Printf("Number of expanded seeds is: %d", len(hits))
	for i := 0; i < len(hits); i++ {
		log.Printf("Full path is : %v\n", sumLen(hits[i]))
		log.Printf("Full path is : %v\n", getSeedPath(hits[i]))
		log.Printf("CheckTail : QueryStart=%d, QueryEnd=%d\n", toTail(hits[i]).QueryStart, toTail(hits[i]).QueryStart+toTail(hits[i]).Length)

	}
	printSeedDevInfo(hits)
	seedHits := getSeqTraversal(gg.Nodes[0], []dna.Base{}, 0, 32)
	for j := 0; j < len(seedHits); j++ {
		fmt.Printf("%s\n", dna.BasesToString(seedHits[j]))
	}
	//seq1 := dna.StringToBases("ATATTAGTGTCTAAAAGCTGTGAAAACACAAACCTACTTTAGTGTAACATTAGTTTTTTTTAAATTTGTAAATCCTTCTGTATGTTTGGATAGTTTTTATACACTTTTCTTTTATTTCAGTTCACCCAAAAGCTTTAAAGGTTGTCTTTAG")
	//seq2 := dna.StringToBases("TCCACCTAGAGGAGCCTGTTCTAGAACCGATAACCACCGTTCAACCTTACCTCCCCTTGTTAATACCGCCTATATACCACCGTCGTCAGCTTACCCTGTGAGGGACTAATAGTAAGCTAAACT-------------------GGTCGAGGTGTAGCGAATGAGGAGGGAAGAAATGGGCTACATTTGGCTAC-AGCAGCGAACACGAATGATGTCTTGAAACGTACATCTGAAGGAGGATTTAGCAGTAAGTACAAAATAGAGTGT")
	//log.Printf("Length of 2 sequences are: %d\n", len(seq1))
}
