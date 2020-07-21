package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
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
	if answer != "g.2A>C|MODERATE|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro" {
		log.Println("Output: ", answer)
		log.Println("Expected: g.2A>C|MODERATE|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro")
		t.Errorf("ERROR: Problem annotating variant")
	}
}

func TestVariantToAnnotationLarge(t *testing.T) {
	testGtf := Read("testdata/KRIT1.gtf")
	testVcf := vcf.Read("testdata/KRIT1.vcf")
	testFasta := fasta.FastaMap(fasta.Read("testdata/hg38_chr7.fa"))
	tree := GenesToIntervalTree(testGtf)
	var err error
	var variant *Variant
	for _, val := range testVcf {
		variant, err = VcfToVariant(val, tree, testFasta)
		if err != nil {
			log.Println(err)
		}
		fmt.Printf("\nANSWER: POS=%d VAR=%s\n", variant.Pos, variant.Info)
		fmt.Printf("OUTPUT: %s\n", VariantToAnnotation(variant, testFasta))
	}
}
