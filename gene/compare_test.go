package gene

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var isEqualTest Gene = Gene{
	id:           "test_gene_id",
	startPos:     0,
	posStrand:    true,
	cdsStarts:    []int{2, 7, 11},
	cdsEnds:      []int{4, 9, 13},
	genomeSeq:    dna.StringToBases("ACATGCACCGTTAACG"),
	cdnaSeq:      dna.StringToBases("ACATGCCGTAACG"),
	featureArray: []Feature{-5, -5, 0, 1, 2, -1, -1, 3, 4, 5, -1, 6, 7, 8, -3, -3},
	orig:         goGeneBackup{},
	utrFive:      subSeq{start: 0, end: 2, seq: dna.StringToBases("AC")},
	utrThree:     subSeq{start: 11, end: 13, seq: dna.StringToBases("CG")},
	codingSeq:    subSeq{start: 2, end: 11, seq: dna.StringToBases("ATGCCGTAA")},
}

var isEqualSubseq = subSeq{
	start: 0,
	end:   1,
	seq:   dna.StringToBases("A"),
}

var empty EffectPrediction = EffectPrediction{
	Consequence: Silent,
	CdnaPos:     0,
	CdnaDist:    0,
	AaRef:       []dna.AminoAcid{},
	AaAlt:       []dna.AminoAcid{},
	StopDist:    0,
}

func TestCompareEqual(t *testing.T) {
	if ok, err := equal(&isEqualTest, &isEqualTest); !ok {
		t.Errorf("ERROR: Equal comparison function should return true when provided the same two Gene structs...\n%s", err)
	}
	if !equalSubSeq(isEqualSubseq, isEqualSubseq) {
		t.Errorf("ERROR: Equal subSeq comparison should return true for two inputs that are exactly the same...")
	}
	if ok, err := equalPred(&empty, &empty); !ok {
		t.Errorf("ERROR: Equal Pred comparison should return true for two inputs empty inputs...\n%s", err)
	}
}
