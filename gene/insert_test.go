package gene

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
)

func TestInsertion(t *testing.T) {
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
		cdsStarts:    []int{2, 7, 13},
		cdsEnds:      []int{4, 11, 15},
		genomeSeq:    dna.StringToBases("ACATGCACCATGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCATGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, 6, 7, -1, 8, 9, 10, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 13, end: 15, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 13, seq: dna.StringToBases("ATGCCATGTAA")},
	}

	_, err = Insertion(gene, 8, []dna.Base{dna.A, dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(gene, &correctPos1); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos2 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 9, 13},
		cdsEnds:      []int{4, 11, 15},
		genomeSeq:    dna.StringToBases("ACATGCATACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Insertion(gene, 5, []dna.Base{dna.A, dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(gene, &correctPos2); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos3 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 9, 13},
		cdsEnds:      []int{4, 11, 15},
		genomeSeq:    dna.StringToBases("ACATGCAATCCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Insertion(gene, 6, []dna.Base{dna.A, dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	if ok, diff := equal(gene, &correctPos3); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos4 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{4, 9, 13},
		cdsEnds:      []int{6, 11, 15},
		genomeSeq:    dna.StringToBases("AATCATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("AATCATGCCGTAACG"),
		featureArray: []Feature{-5, -5, -5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 4, seq: dna.StringToBases("AATC")},
		utrThree:     subSeq{start: 13, end: 15, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 4, end: 13, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Insertion(gene, 0, []dna.Base{dna.A, dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(gene, &correctPos4); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
	Reset(gene)

	var correctPos5 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACATG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTAACATG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3, -3, -3},
		orig:         correctBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 15, seq: dna.StringToBases("CATG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
	}

	_, err = Insertion(gene, 14, []dna.Base{dna.A, dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	if ok, diff := equal(gene, &correctPos5); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
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
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 15},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTCTAACG"),
		cdnaSeq:      dna.StringToBases("ACATGCCGTCTAACG"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, 9, 10, -3, -3},
		orig:         correctNegBackup,
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 13, end: 15, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 13, seq: dna.StringToBases("ATGCCGTCTAA")},
	}

	_, err = Insertion(negGene, 3, []dna.Base{dna.A, dna.G})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}

	if ok, diff := equal(negGene, &correctNeg1); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
}

func TestUndoInsertion(t *testing.T) {
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
		utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
		utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
		codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
		orig:         correctBackup,
	}

	_, _ = Insertion(answerPos, 9, []dna.Base{dna.T})
	Reset(answerPos)
	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing insertion")
	}

	for i := 0; i < len(answerPos.featureArray); i++ {
		if answerPos.featureArray[i] != correctPos.featureArray[i] {
			t.Error("ERROR: Trouble undoing insertion")
		}
	}

	_, _ = Insertion(answerPos, 9, []dna.Base{dna.A, dna.C, dna.T, dna.G})
	_, _ = Insertion(answerPos, 2, []dna.Base{dna.C, dna.T, dna.G})
	_, _ = Insertion(answerPos, 4, []dna.Base{dna.A, dna.C, dna.T})
	_, _ = Insertion(answerPos, 9, []dna.Base{dna.A, dna.C})
	_, _ = Insertion(answerPos, 8, []dna.Base{dna.A, dna.C})
	_, _ = Insertion(answerPos, 0, []dna.Base{dna.A, dna.C})

	Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing multiple insertions")
	}

	for i := 0; i < len(answerPos.featureArray); i++ {
		if answerPos.featureArray[i] != correctPos.featureArray[i] {
			t.Error("ERROR: Trouble undoing multiple insertions")
		}
	}
}

// TODO WIP.
func TestInsertionEffectPrediction(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error
	var pred, correctPred EffectPrediction

	gene := GtfToGene(g["test_gene_id"], f)

	// TEST 1: cDNA Duplication
	_, _ = Insertion(gene, 14, []dna.Base{dna.A, dna.A, dna.A, dna.T, dna.A, dna.T, dna.A, dna.T, dna.A, dna.A, dna.A, dna.A, dna.T})
	pred, err = Insertion(gene, 2, []dna.Base{dna.T, dna.G, dna.C, dna.C})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = Frameshift
	correctPred.CdnaPos = 0
	correctPred.CdnaDist = 0
	correctPred.AaPos = 2
	correctPred.AaRef = []dna.AminoAcid{dna.Stop}
	correctPred.AaAlt = []dna.AminoAcid{dna.Ala}
	correctPred.StopDist = 5

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)

	// TEST 2: Intronic Insertion
	pred, err = Insertion(gene, 5, []dna.Base{dna.T})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = Splice
	correctPred.CdnaPos = 2
	correctPred.CdnaDist = 2
	correctPred.AaPos = 0
	correctPred.AaRef = nil
	correctPred.AaAlt = nil
	correctPred.StopDist = -1

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)

	// TEST 3: Single Base Insertion Causing Frameshift
	pred, err = Insertion(gene, 7, []dna.Base{dna.A})
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
	correctPred.AaAlt = []dna.AminoAcid{dna.His}
	correctPred.StopDist = -2

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)

	// TEST 4: In-frame Insertion
	pred, err = Insertion(gene, 3, []dna.Base{dna.A, dna.A, dna.A})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = InFrameInsertion
	correctPred.CdnaPos = 1
	correctPred.CdnaDist = 0
	correctPred.AaPos = 0
	correctPred.AaRef = []dna.AminoAcid{dna.Met}
	correctPred.AaAlt = []dna.AminoAcid{dna.Ile, dna.Lys}
	correctPred.StopDist = -1

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)

	// TEST 5: Frameshift with stop across exons
	_, _ = Insertion(gene, 7, []dna.Base{dna.T, dna.A, dna.A})
	pred, err = Insertion(gene, 2, []dna.Base{dna.C, dna.C})
	if err != nil {
		if err != ErrNoStopFound {
			t.Error(err)
		}
	}
	correctPred.Consequence = Frameshift
	correctPred.CdnaPos = 0
	correctPred.CdnaDist = 0
	correctPred.AaPos = 0
	correctPred.AaRef = []dna.AminoAcid{dna.Met}
	correctPred.AaAlt = []dna.AminoAcid{dna.Thr}
	correctPred.StopDist = 2

	if ok, diff := equalPred(&pred, &correctPred); !ok {
		t.Errorf("trouble with insertion EffectPrediction. Error is in %s", diff)
	}
	Reset(gene)
}
