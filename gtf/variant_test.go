package gtf

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
)

func TestVcfToVariant(t *testing.T) {
	var inVcf vcf.Vcf
	var inGtf map[string]*Gene = make(map[string]*Gene)
	var inSeq map[string][]dna.Base = make(map[string][]dna.Base)

	inVcf.Chr = "DummyChr"
	inVcf.Ref = "A"
	inVcf.Alt = []string{"C"}
	inVcf.Pos = 2

	inSeq["DummyChr"] = []dna.Base{dna.C, dna.A, dna.T}

	inGtf["DummyGene"] = &Gene{
		GeneID:   "DummyGene",
		GeneName: "DummyGene",
		Transcripts: []*Transcript{{
			Chr:          "DummyChr",
			Strand:       true,
			Start:        1,
			End:          3,
			TranscriptID: "DummyTranscriptID",
			Exons: []*Exon{{
				Start: 1,
				End:   3,
				Cds: &Cds{
					Start: 1,
					End:   3,
					Frame: 0,
				},
			}},
		}},
	}

	// Seq: CAT
	// Variant c.2A>C p.H1P
	tree := GenesToIntervalTree(inGtf)
	variant, err := VcfToVariant(inVcf, tree, inSeq, false)

	if err != nil {
		t.Errorf("%s", err)
	}

	if variant.AaPos != 1 ||
		variant.AaRef[0] != dna.His ||
		variant.AaAlt[0] != dna.Pro ||
		variant.CdnaPos != 2 {
		t.Errorf("ERROR: Problem converting vcf to variant")
	}

}
