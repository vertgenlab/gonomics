package gene

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
)

func TestGtfToGoGene(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var ok bool
	var err error

	// Positive posStrand test
	answerPos := GtfToGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	var correctPos Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	if ok, err = equal(answerPos, &correctPos); !ok {
		t.Errorf("ERROR: Trouble converting gtf to Gene on positive strand, error is %s", err)
	}

	// Negative posStrand test
	answerNeg := GtfToGene(g["test_gene_id_negative"], f)

	var correctNeg Gene = Gene{
		id:           "test_gene_id_negative",
		startPos:     15,
		posStrand:    false,
		cdsStarts:    correctPos.cdsStarts,
		cdsEnds:      correctPos.cdsEnds,
		genomeSeq:    correctPos.genomeSeq,
		cdnaSeq:      correctPos.cdnaSeq,
		featureArray: correctPos.featureArray,
		orig:         correctBackup,
		utrFive:      correctPos.utrFive,
		utrThree:     correctPos.utrThree,
		codingSeq:    correctPos.codingSeq,
	}

	correctNeg.orig.startPos = 15
	if ok, err = equal(answerNeg, &correctNeg); !ok {
		t.Errorf("ERROR: Trouble converting gtf to Gene on negative strand, error is %s", err)
	}
}

func TestPositionConversion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive Tests
	genePos := GtfToGene(g["test_gene_id"], f)
	var testsPass bool = true

	var answerPos, answerOffset int
	var err error

	answerPos, answerOffset, err = GenomicPosToCdna(genePos, 0)
	if answerPos != 0 || answerOffset != -2 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(genePos, 3)
	if answerPos != 1 || answerOffset != 0 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(genePos, 10)
	if answerPos != 5 || answerOffset != 1 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(genePos, 14)
	if answerPos != 8 || answerOffset != 1 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(genePos, 0)
	if answerPos != 2 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(genePos, 3)
	if answerPos != 7 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(genePos, 6)
	if answerPos != 11 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(genePos, 8)
	if answerPos != 13 || err != nil {
		testsPass = false
	}

	if !testsPass {
		t.Error("ERROR: Trouble converting positions on positive posStrand")
	}

	// Negative Tests
	testsPass = true
	geneNeg := GtfToGene(g["test_gene_id_negative"], f)

	answerPos, answerOffset, err = GenomicPosToCdna(geneNeg, 0)
	if answerPos != 8 || answerOffset != 2 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(geneNeg, 3)
	if answerPos != 7 || answerOffset != 0 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(geneNeg, 10)
	if answerPos != 2 || answerOffset != 1 || err != nil {
		testsPass = false
	}

	answerPos, answerOffset, err = GenomicPosToCdna(geneNeg, 14)
	if answerPos != 0 || answerOffset != -1 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(geneNeg, 0)
	if answerPos != 13 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(geneNeg, 3)
	if answerPos != 8 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(geneNeg, 6)
	if answerPos != 4 || err != nil {
		testsPass = false
	}

	answerPos, err = CodingPosToGenomic(geneNeg, 8)
	if answerPos != 2 || err != nil {
		testsPass = false
	}

	if !testsPass {
		t.Error("ERROR: Trouble converting positions on negative posStrand")
	}
}
