package gene

import (
	"errors"
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

	// Positive posStrand test
	posGene := GtfToGoGene(g["test_gene_id"], f)

	answerPos, err = PointMutation(posGene, 6, dna.T)
	if err != nil {
		t.Error(err)
	}

	if answerPos.CdnaPos != 3 ||
		answerPos.CdnaOffset != -1 ||
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
		answerPos.CdnaOffset != 0 ||
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
		answerPos.CdnaOffset != 0 ||
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
		answerPos.CdnaOffset != 0 ||
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
	negGene := GtfToGoGene(g["test_gene_id_negative"], f)

	answerNeg, err = PointMutation(negGene, 9, dna.A)
	if err != nil {
		t.Error(err)
	}

	if answerNeg.CdnaPos != 3 ||
		answerNeg.CdnaOffset != -1 ||
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
		answerNeg.CdnaOffset != 0 ||
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
		answerNeg.CdnaOffset != 0 ||
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
		answerNeg.CdnaOffset != 0 ||
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
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive posStrand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
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

	correctNegBackup := goGeneBackup{
		startPos:     15,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctNeg Gene = Gene{
		id:           "test_gene_id_negative",
		startPos:     15,
		posStrand:    false,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctNegBackup,
	}

	// Negative posStrand test
	answerNeg := GtfToGoGene(g["test_gene_id_negative"], f)

	_, _ = PointMutation(answerNeg, 9, dna.G)
	Reset(answerNeg)

	if ok, diff := equal(&correctNeg, answerNeg); !ok {
		t.Errorf("trouble undoing point mutation on negative posStrand. Problem with %s", diff)
	}
}

func TestUndoInsertion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive posStrand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
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

func TestUndoDeletion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")

	// Positive posStrand test
	answerPos := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
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

	_, _ = Deletion(answerPos, 9, 13)
	_, _ = Deletion(answerPos, 2, 5)
	_, _ = Deletion(answerPos, 4, 6)

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

func TestInsertion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error

	gene := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos1 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 7, 13},
		cdsEnds:      []int{4, 11, 15},
		genomeSeq:    dna.StringToBases("ACATGCACCATGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCATGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, 6, 7, -1, 8, 9, 10, -3, -3},
		orig:         correctBackup,
	}

	_, err = Insertion(gene, 8, []dna.Base{dna.A, dna.T})
	if err != nil {
		t.Error(err)
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
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
	}

	_, err = Insertion(gene, 5, []dna.Base{dna.A, dna.T})
	if err != nil {
		t.Error(err)
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
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
		orig:         correctBackup,
	}

	_, err = Insertion(gene, 6, []dna.Base{dna.A, dna.T})
	if err != nil {
		t.Error(err)
	}
	if ok, diff := equal(gene, &correctPos3); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}

	// Negative posStrand test
	negGene := GtfToGoGene(g["test_gene_id_negative"], f)

	correctNegBackup := goGeneBackup{
		startPos:     15,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctNeg1 Gene = Gene{
		id:           "test_gene_id_negative",
		startPos:     15,
		posStrand:    false,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 15},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTCTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTCTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, 9, 10, -3, -3},
		orig:         correctNegBackup,
	}

	_, err = Insertion(negGene, 3, []dna.Base{dna.A, dna.G})
	if err != nil {
		t.Error(err)
	}

	if ok, diff := equal(negGene, &correctNeg1); !ok {
		t.Errorf("trouble making insertion. Error is in %s", diff)
	}
}

func TestDeletion(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error

	gene := GtfToGoGene(g["test_gene_id"], f)

	correctBackup := goGeneBackup{
		startPos:     0,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctPos1 Gene = Gene{
		id:           "test_gene_id",
		startPos:     0,
		posStrand:    true,
		cdsStarts:    []int{2, 4},
		cdsEnds:      []int{3, 5},
		genomeSeq:    dna.StringToBases("ACATAACG"),
		cdnaSeq:      dna.StringToBases("ATAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, 3, -3, -3},
		orig:         correctBackup,
	}

	_, err = Deletion(gene, 3, 12)
	if err != nil {
		t.Error(err)
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
		cdnaSeq:      dna.StringToBases("CCGTAA"),
		featureArray: []Feature{-5, -5, -1, -1, 0, 1, 2, -1, 3, 4, 5, -3, -3},
		orig:         correctBackup,
	}

	_, err = Deletion(gene, 1, 5)
	if err != nil {
		t.Error(err)
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
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, 6, 7, 8, -3, -3},
		orig:         correctBackup,
	}

	_, err = Deletion(gene, 9, 11)
	if err != nil {
		t.Error(err)
	}
	if ok, diff := equal(gene, &correctPos3); !ok {
		t.Errorf("trouble making deletion. Error is in %s", diff)
	}

	// Negative posStrand test
	negGene := GtfToGoGene(g["test_gene_id_negative"], f)

	correctNegBackup := goGeneBackup{
		startPos:     15,
		cdsStarts:    []int{2, 7, 11},
		cdsEnds:      []int{4, 9, 13},
		genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
		cdnaSeq:      dna.StringToBases("ATGCCGTAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	}

	var correctNeg1 Gene = Gene{
		id:           "test_gene_id_negative",
		startPos:     15,
		posStrand:    false,
		cdsStarts:    []int{2, 4},
		cdsEnds:      []int{3, 5},
		genomeSeq:    dna.StringToBases("ACATAACG"),
		cdnaSeq:      dna.StringToBases("ATAA"),
		featureArray: []Feature{-5, -5, 0, 1, 2, 3, -3, -3},
		orig:         correctNegBackup,
	}

	_, err = Deletion(negGene, 3, 12)
	if err != nil {
		t.Error(err)
	}

	if ok, diff := equal(negGene, &correctNeg1); !ok {
		t.Errorf("trouble making deletion on negative posStrand. Error is in %s", diff)
	}
}

func equal(alpha, beta *Gene) (bool, error) {
	if alpha.id != beta.id {
		return false, errors.New("id")
	}

	if alpha.startPos != beta.startPos {
		return false, errors.New("startPos")
	}

	if alpha.posStrand != beta.posStrand {
		return false, errors.New("posStrand")
	}

	if len(alpha.cdsStarts) != len(beta.cdsStarts) {
		return false, errors.New("cdsStarts")
	}

	for idx, val := range alpha.cdsStarts {
		if val != beta.cdsStarts[idx] {
			return false, errors.New("cdsStarts")
		}
	}

	if len(alpha.cdsEnds) != len(beta.cdsEnds) {
		return false, errors.New("cdsEnds")
	}

	for idx, val := range alpha.cdsEnds {
		if val != beta.cdsEnds[idx] {
			return false, errors.New("cdsEnds")
		}
	}

	if len(alpha.genomeSeq) != len(beta.genomeSeq) {
		return false, errors.New("genomeSeq")
	}

	for idx, val := range alpha.genomeSeq {
		if val != beta.genomeSeq[idx] {
			return false, errors.New("genomeSeq")
		}
	}

	if len(alpha.cdnaSeq) != len(beta.cdnaSeq) {
		return false, errors.New("cdnaSeq")
	}

	for idx, val := range alpha.cdnaSeq {
		if val != beta.cdnaSeq[idx] {
			return false, errors.New("cdnaSeq")
		}
	}

	if len(alpha.featureArray) != len(beta.featureArray) {
		return false, errors.New("featureArray")
	}

	for idx, val := range alpha.featureArray {
		if val != beta.featureArray[idx] {
			return false, errors.New("featureArray")
		}
	}

	return true, nil
}

func TestInsertionEffectPrediction(t *testing.T) {
	g := gtf.Read("testdata/test.gtf")
	f := fasta.Read("testdata/test.fasta")
	var err error
	var pred EffectPrediction

	gene := GtfToGoGene(g["test_gene_id"], f)

	_, _ = Insertion(gene, 14, []dna.Base{dna.A, dna.A, dna.A, dna.T, dna.A, dna.T, dna.A, dna.A, dna.A, dna.T, dna.A, dna.A, dna.T})

	pred, err = Insertion(gene, 2, []dna.Base{dna.A})
	if err != nil {
		t.Error(err)
	}

	printEffPred(pred)
}

func printEffPred(pred EffectPrediction) {
	var consequence string
	switch pred.Consequence {
	case Intronic:
		consequence = "Intronic"
	case Silent:
		consequence = "Silent"
	case Missense:
		consequence = "Missense"
	case Nonsense:
		consequence = "Nonsense"
	case Frameshift:
		consequence = "Frameshift"
	case Intergenic:
		consequence = "Intergenic"
	case Splice:
		consequence = "Splice"
	case FarSplice:
		consequence = "FarSplice"
	case DisruptStart:
		consequence = "DisruptStart"
	case DisruptStop:
		consequence = "DisruptStop"
	}
	fmt.Printf("Consequence: %s\n", consequence)
	fmt.Printf("cDNA Pos: %d%+d\n", pred.CdnaPos, pred.CdnaOffset)
	fmt.Printf("AAPos: %d\n", pred.AaPos)
	fmt.Printf("AaRef: %s\n", dna.PolypeptideToString(pred.AaRef))
	fmt.Printf("AaAlt: %s\n", dna.PolypeptideToString(pred.AaAlt))
}

//TODO Indel EffectPred tests
