package gene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
	"testing"
)

func TestDeletion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error

	gene := GtfToGene(g["test_gene_id"], f)

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

	var correctPos1 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 4},
		cdsEnds:      []int{3, 5},
		genomeSeq:    dna.StringToBases("ACATAACG"),
		cdnaSeq:      dna.StringToBases("ACATAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, 3, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 6, end: 8, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 6, seq: dna.StringToBases("ATAA")},
	}

	_, err = Deletion(gene, 4, 12)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(gene, &correctPos1); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos2 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{4, 8},
		cdsEnds:      []int{6, 10},
		genomeSeq:    dna.StringToBases("ACCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACCCGTAACG"),
		featureArray: []Feature{-5, -5, -1, -1, 0, 1, 2, -1, 3, 4, 5, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 8, end: 10, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 8, seq: dna.StringToBases("CCGTAA")},
	}

	_, err = Deletion(gene, 2, 5)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(gene, &correctPos2); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos3 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 10},
		cdsEnds:      []int{4, 9, 12},
		genomeSeq:    dna.StringToBases("ACATGCACCGTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Deletion(gene, 10, 11)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	if ok, diff := equal(gene, &correctPos3); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos4 Gene = Gene{
		id:           "test_gene_id",
		startPos:     2,
		posStrand:    true,
		cdsStarts:    []int{0, 5, 9},
		cdsEnds:      []int{2, 7, 11},
		genomeSeq:    dna.StringToBases("ATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAACG"),
		featureArray: []Feature{0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 0, seq: dna.StringToBases("")},
		utrThree:     subSeq{start: 9, end: 11, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 0, end: 9, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Deletion(gene, 0, 2)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	if ok, diff := equal(gene, &correctPos4); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}

	Reset(gene)

	var correctPos5 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAA"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 11, seq: dna.StringToBases("")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Deletion(gene, 14, 16)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	if ok, diff := equal(gene, &correctPos5); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}

	// Negative posStrand test
	negGene := GtfToGene(g["test_gene_id_negative"], f)

	correctNegBackup := goGeneBackup{
		startPos:     15,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	var correctNeg1 Gene = Gene{
		id:           "test_gene_id_negative",
		startPos:     15,
		posStrand:    false,
		cdsStarts:    []int{2, 4},
		cdsEnds:      []int{3, 5},
		genomeSeq:    dna.StringToBases("ACATAACG"),
		cdnaSeq:      dna.StringToBases("ACATAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, 3, -3, -3},
		orig:         correctNegBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 6, end: 8, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 6, seq: dna.StringToBases("ATAA")},
	}

	_, err = Deletion(negGene, 4, 12)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(negGene, &correctNeg1); !ok {
		t.Errorf("trouble making deletion on negative posStrand. Error is in %s", diff)
	}
}

func TestUndoDeletion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

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

	_, _ = Deletion(answerPos, 10, 13)
	Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing deletion")
	}

	for i := 0; i < len(answerPos.featureArray); i++ {
		if answerPos.featureArray[i] != correctPos.featureArray[i] {
			t.Error("ERROR: Trouble undoing deletion")
		}
	}

	_, _ = Deletion(answerPos, 10, 13)
	_, _ = Deletion(answerPos, 3, 5)
	_, _ = Deletion(answerPos, 5, 6)

	Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing multiple deletions")
	}

	for i := 0; i < len(answerPos.featureArray); i++ {
		if answerPos.featureArray[i] != correctPos.featureArray[i] {
			t.Error("ERROR: Trouble undoing multiple deletions")
		}
	}
}

func TestDeletionEffectPrediction(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error
	var pred, correctPred EffectPrediction

	gene := GtfToGene(g["test_gene_id"], f)

	// TEST 1: In-Frame Exon Deletion
	pred, err = Deletion(gene, 7, 10)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = InFrameDeletion
	correctPred.CdnaPos = 3
	correctPred.CdnaDist = 0
	correctPred.AaPos = 1
	correctPred.AaRef = []dna.AminoAcid{dna.Pro}
	correctPred.AaAlt = []dna.AminoAcid{}
	correctPred.StopDist = -1

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)

	// TEST 1: 1bp Deletion in CDS
	pred, err = Deletion(gene, 7, 8)
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = Frameshift
	correctPred.CdnaPos = 3
	correctPred.CdnaDist = 0
	correctPred.AaPos = 1
	correctPred.AaRef = []dna.AminoAcid{dna.Pro}
	correctPred.AaAlt = []dna.AminoAcid{}
	correctPred.StopDist = -1

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)
}
