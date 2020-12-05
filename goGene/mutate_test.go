package goGene

import (
	"fmt"
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
	Reset(posGene)
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
	Reset(posGene)
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
	Reset(posGene)
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
	Reset(posGene)
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
	Reset(negGene)
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
	Reset(negGene)
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
	Reset(negGene)
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
	Reset(negGene)
	if err != nil {
		t.Error(err)
	}
}

func TestUndoPointMutation(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive strand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos GoGene = GoGene{
		id:           "test_gene_id",
		startPos:     0,
		strand:       true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
	}

	_, _ = PointMutation(answerPos, 9, dna.T)
	Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing point mutation")
	}

	_, _ = PointMutation(answerPos, 9, dna.T)
	_, _ = PointMutation(answerPos, 2, dna.C)
	_, _ = PointMutation(answerPos, 4, dna.A)
	_, _ = PointMutation(answerPos, 9, dna.G)
	Reset(answerPos)

	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 {
		t.Error("ERROR: Trouble undoing multiple point mutations")
	}
}

//func TestUndoInsertion(t *testing.T) {
//	g := gtf.Read("testdata/test.gtf")
//	f := fasta.Read("testdata/test.fasta")
//
//	// Positive strand test
//	answerPos := GtfToGoGene(g["test_gene_id"], f)
//
//	var correctPos GoGene = GoGene{
//		id:           "test_gene_id",
//		startPos:     0,
//		strand:       true,
//		cdsStarts:    []int{2, 7, 11},
//		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
//		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
//		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
//      origFeatureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
//	}
//
//	_, _ = Insertion(answerPos, 9, []dna.Base{dna.T})
//	_, err := Undo(answerPos)
//	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
//		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 || err != nil {
//		t.Error("ERROR: Trouble undoing insertion")
//	}
//
//	for i := 0; i < len(answerPos.featureArray); i++ {
//		if answerPos.featureArray[i] != correctPos.featureArray[i] {
//			t.Error("ERROR: Trouble undoing insertion")
//		}
//	}
//
//	_, _ = Insertion(answerPos, 9, []dna.Base{dna.A, dna.C, dna.T, dna.G})
//	_, _ = Insertion(answerPos, 2, []dna.Base{dna.C, dna.T, dna.G})
//	_, _ = Insertion(answerPos, 4, []dna.Base{dna.A, dna.C, dna.T})
//	_, _ = Insertion(answerPos, 9, []dna.Base{dna.A, dna.C})
//	_, _ = Insertion(answerPos, 8, []dna.Base{dna.A, dna.C})
//	_, _ = Insertion(answerPos, 0, []dna.Base{dna.A, dna.C})
//
//	numUndone, err := Reset(answerPos)
//
//	if dna.CompareSeqsIgnoreCase(correctPos.genomeSeq, answerPos.genomeSeq) != 0 ||
//		dna.CompareSeqsIgnoreCase(correctPos.cdnaSeq, answerPos.cdnaSeq) != 0 ||
//		err != nil || numUndone != 6 {
//		t.Error("ERROR: Trouble undoing multiple insertions")
//	}
//
//	for i := 0; i < len(answerPos.featureArray); i++ {
//		if answerPos.featureArray[i] != correctPos.featureArray[i] {
//			t.Error("ERROR: Trouble undoing multiple insertions")
//		}
//	}
//}

func TestUndoDeletion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive strand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos GoGene = GoGene{
		id:           "test_gene_id",
		startPos:     0,
		strand:       true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
	}
	_, _ = Deletion(answerPos, 9, 13)
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

	fmt.Println(answerPos)
	_, _ = Deletion(answerPos, 9, 13)
	fmt.Println(answerPos)
	_, _ = Deletion(answerPos, 2, 5)
	fmt.Println(answerPos)
	_, _ = Deletion(answerPos, 4, 6)
	fmt.Println(answerPos)

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

//func TestInsertion(t *testing.T) {
//	g := gtf.Read("testdata/test.gtf")
//	f := fasta.Read("testdata/test.fasta")
//
//	// Positive strand test
//	posGene := GtfToGoGene(g["test_gene_id"], f)
//	fmt.Println(posGene)
//
//	//fmt.Println(dna.BasesToString(posGene.genomeSeq))
//	//fmt.Println(posGene.featureArray)
//	//fmt.Println(dna.BasesToString(posGene.cdnaSeq))
//	Insertion(posGene, 8, []dna.Base{dna.A, dna.T})
//	//fmt.Println(dna.BasesToString(posGene.genomeSeq))
//	//fmt.Println(posGene.featureArray)
//	//fmt.Println(dna.BasesToString(posGene.cdnaSeq))
//
//	fmt.Println(posGene)
//}

//func TestDeletion(t *testing.T) {
//	g := gtf.Read("testdata/test.gtf")
//	f := fasta.Read("testdata/test.fasta")
//
//	// Positive strand test
//	posGene := GtfToGoGene(g["test_gene_id"], f)
//
//	fmt.Println(posGene)
//
//	fmt.Println(dna.BasesToString(posGene.genomeSeq))
//	fmt.Println(posGene.featureArray)
//	fmt.Println(dna.BasesToString(posGene.cdnaSeq))
//	Deletion(posGene, 3, 12)
//	fmt.Println(dna.BasesToString(posGene.genomeSeq))
//	fmt.Println(posGene.featureArray)
//	fmt.Println(dna.BasesToString(posGene.cdnaSeq))
//
//	fmt.Println(posGene)
//}
