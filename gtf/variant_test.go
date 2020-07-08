package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
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
	answer, err := VcfToVariant(inVcf, tree, inSeq)
	common.ExitIfError(err)

	fmt.Println(answer)
}
