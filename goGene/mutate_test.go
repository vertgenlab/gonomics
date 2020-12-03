package goGene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
	"testing"
)

func TestPointMutation(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	var answerPos, answerNeg EffectPrediction
	var err error

	// Positive strand test
	posGene := GtfToGoGene(g["test_gene_id"], f)

	answerPos, err = PointMutation(posGene, 6, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 3 ||
		answerPos.CdnaOffset != -1 ||
		answerPos.Consequence != Splice {
		t.Error("trouble with intronic point mutation on positive strand")
	}
	_, err = Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 3, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 1 ||
		answerPos.CdnaOffset != 0 ||
		answerPos.Consequence != DisruptStart ||
		answerPos.AaPos != 0 ||
		answerPos.AaRef[0] != dna.Met ||
		answerPos.AaAlt[0] != dna.Lys {
		t.Error("trouble with start point mutation on positive strand")
	}
	_, err = Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 8, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 4 ||
		answerPos.CdnaOffset != 0 ||
		answerPos.Consequence != Missense ||
		answerPos.AaPos != 1 ||
		answerPos.AaRef[0] != dna.Pro ||
		answerPos.AaAlt[0] != dna.Gln {
		t.Error("trouble with missense point mutation on positive strand")
	}
	_, err = Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 11, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 6 ||
		answerPos.CdnaOffset != 0 ||
		answerPos.Consequence != DisruptStop ||
		answerPos.AaPos != 2 ||
		answerPos.AaRef[0] != dna.Stop ||
		answerPos.AaAlt[0] != dna.Lys {
		t.Error("trouble with disrupt stop point mutation on positive strand")
	}
	_, err = Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	// Negative strand test
	negGene := GtfToGoGene(g["test_gene_id_negative"], f)

	answerNeg, err = PointMutation(negGene, 9, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 3 ||
		answerNeg.CdnaOffset != -1 ||
		answerNeg.Consequence != Splice {
		t.Error("trouble with intronic point mutation on negative strand")
	}
	_, err = Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 12, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 1 ||
		answerNeg.CdnaOffset != 0 ||
		answerNeg.Consequence != DisruptStart ||
		answerNeg.AaPos != 0 ||
		answerNeg.AaRef[0] != dna.Met ||
		answerNeg.AaAlt[0] != dna.Lys {
		t.Error("trouble with start point mutation on negative strand")
	}
	_, err = Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 7, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 4 ||
		answerNeg.CdnaOffset != 0 ||
		answerNeg.Consequence != Missense ||
		answerNeg.AaPos != 1 ||
		answerNeg.AaRef[0] != dna.Pro ||
		answerNeg.AaAlt[0] != dna.Gln {
		t.Error("trouble with missense point mutation on negative strand")
	}
	_, err = Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 4, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 6 ||
		answerNeg.CdnaOffset != 0 ||
		answerNeg.Consequence != DisruptStop ||
		answerNeg.AaPos != 2 ||
		answerNeg.AaRef[0] != dna.Stop ||
		answerNeg.AaAlt[0] != dna.Lys {
		t.Error("trouble with disrupt stop point mutation on negative strand")
	}
	_, err = Reset(negGene)
	if err != nil {
		t.Error(err)
	}
}

func TestUndoPointMutation(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive strand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	var correctPos GoGene = GoGene{
		id:           "test_gene_id",
		startPos:     0,
		strand:       true,
		cdsStarts:    []int{2, 7, 11},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	_, _ = PointMutation(answerPos, 9, dna.T)
	_, err := Undo(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 || err != nil {
		t.Error("ERROR: Trouble undoing point mutation")
	}

	_, _ = PointMutation(answerPos, 9, dna.T)
	_, _ = PointMutation(answerPos, 2, dna.C)
	_, _ = PointMutation(answerPos, 4, dna.A)
	_, _ = PointMutation(answerPos, 9, dna.G)

	numUndone, err := Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 ||
		err != nil || numUndone != 4 {
		t.Error("ERROR: Trouble undoing multiple point mutations")
	}
}
