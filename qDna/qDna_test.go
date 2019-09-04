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

func TestIndex(t *testing.T) {
	fastaFile := FromFastaSlice(fasta.Read("testdata/chrM.fa"))
	indexFile := IndexRefSlidingWindow(fastaFile, 30)
	Write("testdata/chrM_map.ham5", indexFile)
	//Read("testdata/chrM_map.ham5")
}

func TestIndexing(t *testing.T) {
	var bases = dna.StringToBases("ATGATGTATGTATGTATGT")
	var bases2 = dna.StringToBases("ATG")
	loc := Location{Assembly: "", Chr: "chrI", Start: 0, End: 0}
	var testFrags []*QFrag
	testFrag := &QFrag{Seq: FromDna(bases), From: []*Location{&loc}, Fwd: nil,  Rev: nil}
	testFrags = append(testFrags, testFrag)
	loc2 := Location{Assembly: "", Chr: "chrII", Start: 0, End: 0}
	testFrag2 := &QFrag{Seq: FromDna(bases2), From: []*Location{&loc2}, Fwd: nil, Rev: nil}
	testFrags = append(testFrags, testFrag2)

	m := indexRef(testFrags, 3)

	for keys := range m {
		fmt.Println(m[keys])
	}
	fmt.Println("position: ", m[putTogether(dna.StringToBases("ATG"))])
	//alpha := &NoGapAln{Start: 100, End: 200, Score: 0}
	//beta := &NoGapAln{Start: 10, End: 20, Score: 0}
	//merged := mergeTwoSeed(alpha, beta)
	//for i := 0; i < len(merged); i++ {
	//	fmt.Println(merged[i])
	//}
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

func TestAlign(t *testing.T) {
	fastaFile := FromFastaSlice(fasta.Read("testdata/chrM.fa"))
	fastq := fastq.Read("testdata/CL13_chrM.200.fastq")

	fmt.Println("Sam 200: ")
	start1 := time.Now()
	sam := GSW2(fastaFile, fastq)
	t1 := time.Since(start1)

	//for j := 0; j < len(sam); j++ {
	fmt.Println(sam[0])
	//}
	fmt.Println("GSW Aligner v2\nfaster.better,stronger: ", t1)
	fmt.Println("Number of reads: ", len(fastq))
	fmt.Println("Number of sam records: ", len(sam))

}

/*
func TestAlign2(t *testing.T) {
	fastaFile := FromFastaSlice(fasta.Read("testdata/chrM.fa"))
	fastq2 := fastq.Read("testdata/CL13_chrM.mapped.1000.fastq")
	fmt.Println("Sam test Mapped: ")
	start2 := time.Now()
	sam2 := GSW2(fastaFile, fastq2)
	t2 := time.Since(start2)

	//for j := 0; j < len(sam2); j++ {
	fmt.Println(sam2[0])
	//}
	fmt.Println("GSW Aligner v2\nfaster.better,stronger: ", t2)
}

/*func TestQDnaScoreLoop(t *testing.T) {
	fq := fastq.Read("testdata/CL13_chrM.unmapped.100.fastq")
	var testScore float64
	var testScore2 float64
	var testScore3 float64
	start := time.Now()
	for i := 0; i < len(fq); i++ {
		for j := 0; j < len(fq[i].Seq); j++ {
			testScore += QDnaFasterScore(FromFastq(fq[i])[j], FromFastq(fq[i])[j], HumanChimpTwoScoreMatrix)
		}

	}
	t1 := time.Since(start)
	//fmt.Println("Craig's version: ", t1)
	start2 := time.Now()
	for i := 0; i < len(fq); i++ {
		for j := 0; j < len(fq[i].Seq); j++ {
			testScore2 += QDnaScore(FromFastq(fq[i])[j], FromFastq(fq[i])[j], HumanChimpTwoScoreMatrix)
		}
	}
	t2 := time.Since(start2)
}*/

/*
func TestConvert(t *testing.T) {
	for _, test := range convertTests {
		actual := fasta.Read(test.filename)
		qfrag := FromFastaSlice(actual)
		qfragOutput := toFastaList(qfrag)
		if !fasta.AllAreEqual(qfragOutput, actual) {
			t.Errorf("The %s file did not convert properly.", test.filename)
		}

		//a := &QBase{A: 1, C: 0, G: 0, T: 0}
		//b := &QBase{A: 0, C: 1, G: 0, T: 0}
		//fmt.Println(QDnaScore(a, b, HumanChimpTwoScoreMatrix))
		//fmt.Println("DONE")
	}

}*/
