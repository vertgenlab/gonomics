package qDna

import (
	"fmt"
	//"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"testing"
	"time"
)

var convertTests = []struct {
	filename string // input
}{
	{"testdata/seq.fa"},
}

var bases = dna.StringToBases("ATGATGG")
var testFasta = "testdata/seq.fa"
var printTest = QFrag{Seq: FromDna(bases), From: nil, Fwd: nil, Rev: nil}

func TestAlign(t *testing.T) {
	fastaFile := FromFastaSlice(fasta.Read("testdata/chrM.fa"))
	fastq := fastq.Read("testdata/CL13_chrM.200.fastq")
	indexFile := IndexRefSlidingWindow(fastaFile, 20)
	Write("testdata/chrM_map.ham5", indexFile)
	ham5 := Read("testdata/chrM_map.ham5")

	fmt.Println("Sam 200: ")
	start1 := time.Now()
	sam := GSW(fastaFile, fastq, ham5)
	t1 := time.Since(start1)

	//for j := 0; j < len(sam); j++ {
	fmt.Println(sam[0])
	//}
	fmt.Println("GSW Aligner v2\nfaster.better,stronger: ", t1)
	fmt.Println("Number of reads: ", len(fastq))
	fmt.Println("Number of sam records: ", len(sam))

}

/*
func TestPrint(t *testing.T) {
	printFasta := fasta.Read(testFasta)
	//fmt.Println(printFasta)
	seq1 := FromFasta(printFasta[0])
	seq2 := FromFasta(printFasta[1])
	fmt.Println("Sequence 1")
	for i := 0; i < len(seq1.Seq); i++ {
		fmt.Println(seq1.Seq[i])
	}
	fmt.Println("Sequence 2")
	for j := 0; j < len(seq2.Seq); j++ {
		fmt.Println(seq2.Seq[j])
	}
	score, alignment, _, maxI, _, _ := SmithWaterman(seq1.Seq, seq2.Seq, HumanChimpTwoScoreMatrix, -600)
	fmt.Println(score, alignment, maxI)
	fmt.Println(align.LocalView(mostLikelySeq(seq1.Seq), mostLikelySeq(seq2.Seq), alignment, maxI))
	//fmt.Println(LocalMaxScore(-1,-1100,-91))
	//fmt.Println("Testing Affine Gap: \n\n ")
	//blastz, alignment2 := AffineGap(seq1.Seq, seq2.Seq, HumanChimpTwoScoreMatrix, 600, 150)
	//fmt.Println(blastz, alignment2)
	//fmt.Println(align.View(mostLikelySeq(seq1.Seq), mostLikelySeq(seq2.Seq), alignment2))

}*/
