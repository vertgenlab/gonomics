package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
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
			Strand: 	true,
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
	if answer != "g.DummyChr:2A>C|Missense|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro" {
		log.Println("Output: ", answer)
		log.Println("Expected: g.DummyChr:2A>C|Missense|DummyGene|DummyTranscriptID:c.2A>C|p.His1Pro")
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
	var words, newWords, newerWords []string
	var correctString, annotation string
	var errorCount int
	for _, val := range testVcf {
		variant, err = VcfToVariant(val, tree, testFasta)
		if err != nil {
			log.Println(err)
		}
		words = strings.Split(val.Info,"|")
		correctString = words[0]
		annotation = VariantToAnnotation(variant, testFasta)
		newWords = strings.Split(annotation, "|")
		newerWords = strings.Split(newWords[3], ":")
		//fmt.Sprintln(newerWords)
		if newWords[4] == correctString || newerWords[1] == correctString ||
			strings.HasPrefix(correctString, "c.-") || strings.HasPrefix(correctString, "c.*") {
			continue
		}
		errorCount++
		//if correctString == "p.Leu93Pro" { // TODO: REMOVE THIS LINE
			fmt.Printf("\nWARNING: ANNOTATION MISMATCH\n")
			fmt.Printf("EXPECTED: %s\n", correctString)
			fmt.Printf("RECEIVED: %s\n", annotation)
		//}
	}
	if errorCount != 0 {
		t.Errorf("ERROR: %d variants were misannotated", errorCount)
	}
}
