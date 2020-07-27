package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"testing"
)

func TestVcfToVariant(t *testing.T) {
	var inVcf *vcf.Vcf = new(vcf.Vcf)
	var inGtf map[string]*Gene = make(map[string]*Gene)
	var inSeq map[string][]dna.Base = make(map[string][]dna.Base)

	inVcf.Chr = "DummyChr"
	inVcf.Ref = "A"
	inVcf.Alt = "C"
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
				Cds: &CDS{
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
	variant, err := VcfToVariant(inVcf, tree, inSeq)

	if err != nil {
		t.Errorf("%s", err)
	}

	answer := fmt.Sprint(VariantToAnnotation(variant, inSeq))
	if answer != "GoEP=g.DummyChr:2A>C|Missense|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro" {
		log.Println("Output: ", answer)
		log.Println("Expected: g.DummyChr:2A>C|Missense|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro")
		t.Errorf("ERROR: Problem annotating variant")
	}
}
