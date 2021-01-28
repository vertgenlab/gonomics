package gene

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

	// Positive posStrand test
	posGene := GtfToGene(g["test_gene_id"], f)

	answerPos, err = PointMutation(posGene, 6, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 3 ||
		answerPos.CdnaDist != -1 ||
		answerPos.Consequence != Splice {
		t.Error("trouble with intronic point mutation on positive posStrand")
	}
	Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 3, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 1 ||
		answerPos.CdnaDist != 0 ||
		answerPos.Consequence != DisruptStart ||
		answerPos.AaPos != 0 ||
		answerPos.AaRef[0] != dna.Met ||
		answerPos.AaAlt[0] != dna.Lys {
		t.Error("trouble with start point mutation on positive posStrand")
	}
	Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 8, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 4 ||
		answerPos.CdnaDist != 0 ||
		answerPos.Consequence != Missense ||
		answerPos.AaPos != 1 ||
		answerPos.AaRef[0] != dna.Pro ||
		answerPos.AaAlt[0] != dna.Gln {
		t.Error("trouble with missense point mutation on positive posStrand")
	}
	Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	answerPos, err = PointMutation(posGene, 11, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 6 ||
		answerPos.CdnaDist != 0 ||
		answerPos.Consequence != DisruptStop ||
		answerPos.AaPos != 2 ||
		answerPos.AaRef[0] != dna.Stop ||
		answerPos.AaAlt[0] != dna.Lys {
		t.Error("trouble with disrupt stop point mutation on positive posStrand")
	}
	Reset(posGene)
	if err != nil {
		t.Error(err)
	}

	// Negative posStrand test
	negGene := GtfToGene(g["test_gene_id_negative"], f)

	answerNeg, err = PointMutation(negGene, 9, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 3 ||
		answerNeg.CdnaDist != -1 ||
		answerNeg.Consequence != Splice {
		t.Error("trouble with intronic point mutation on negative posStrand")
	}
	Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 12, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 1 ||
		answerNeg.CdnaDist != 0 ||
		answerNeg.Consequence != DisruptStart ||
		answerNeg.AaPos != 0 ||
		answerNeg.AaRef[0] != dna.Met ||
		answerNeg.AaAlt[0] != dna.Lys {
		t.Error("trouble with start point mutation on negative posStrand")
	}
	Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 7, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 4 ||
		answerNeg.CdnaDist != 0 ||
		answerNeg.Consequence != Missense ||
		answerNeg.AaPos != 1 ||
		answerNeg.AaRef[0] != dna.Pro ||
		answerNeg.AaAlt[0] != dna.Gln {
		t.Error("trouble with missense point mutation on negative posStrand")
	}
	Reset(negGene)
	if err != nil {
		t.Error(err)
	}

	answerNeg, err = PointMutation(negGene, 4, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 6 ||
		answerNeg.CdnaDist != 0 ||
		answerNeg.Consequence != DisruptStop ||
		answerNeg.AaPos != 2 ||
		answerNeg.AaRef[0] != dna.Stop ||
		answerNeg.AaAlt[0] != dna.Lys {
		t.Error("trouble with disrupt stop point mutation on negative posStrand")
	}
	Reset(negGene)
	if err != nil {
		t.Error(err)
	}
}

func TestUndoPointMutation(t *testing.T) {
	var ok bool
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

	_, _ = PointMutation(answerPos, 9, dna.T)
	Reset(answerPos)

	if ok, _ = equal(&correctPos, answerPos); !ok {
		t.Error("ERROR: Trouble undoing point mutation")
	}

	_, _ = PointMutation(answerPos, 9, dna.T)
	_, _ = PointMutation(answerPos, 2, dna.C)
	_, _ = PointMutation(answerPos, 4, dna.A)
	_, _ = PointMutation(answerPos, 9, dna.G)
	Reset(answerPos)

	if ok, _ = equal(&correctPos, answerPos); !ok {
		t.Error("ERROR: Trouble undoing multiple point mutations")
	}

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

	// Negative posStrand test
	answerNeg := GtfToGene(g["test_gene_id_negative"], f)

	_, _ = PointMutation(answerNeg, 9, dna.G)
	Reset(answerNeg)

	if ok, diff := equal(&correctNeg, answerNeg); !ok {
		t.Errorf("trouble undoing point mutation on negative posStrand. Problem with %s", diff)
	}
}
