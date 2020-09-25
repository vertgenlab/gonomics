package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
	"testing"
)

func TestVariantToAnnotationLarge(t *testing.T) {
	testGtf := Read("testdata/test.gtf")
	testVcf := vcf.Read("testdata/test.vcf")
	krit := fasta.Read("testdata/krit1.fa")
	cftr := fasta.Read("testdata/cftr.fa")

	// The fasta files were split up to reduce the filesize. The following 6 lines
	// assemble the fasta files so that everything is in the correct place
	f := []*fasta.Fasta{{"chr7", make([]dna.Base, 92198968)}}
	f[0].Seq = append(f[0].Seq, krit[0].Seq...)
	for i := 0; i < 117480024-92246100; i++ {
		f[0].Seq = append(f[0].Seq, dna.N)
	}
	f[0].Seq = append(f[0].Seq, cftr[0].Seq...)

	//f := fasta.Read("testdata/hg38.fa")
	//testGtf := Read("testdata/TTN.gtf")
	//testVcf := vcf.Read("testdata/TTN.vcf")

	fasta.AllToUpper(f)
	testFasta := fasta.FastaMap(f)
	tree := GenesToIntervalTree(testGtf)
	var variant *vcfEffectPrediction
	var words, newWords, newerWords []string
	var correctCDNA, correctProt, outputCDNA, outputProt, annotation string
	var errorCount int
	for _, val := range testVcf {
		variant, _ = VcfToVariant(val, tree, testFasta, false)
		words = strings.Split(val.Info, "|")
		correctCDNA = words[0]
		correctProt = words[1]
		annotation = VariantToAnnotation(variant, testFasta)
		newWords = strings.Split(annotation, "|")
		newerWords = strings.Split(newWords[2], ":")
		outputCDNA = newerWords[1]
		outputProt = newWords[3]
		if (outputCDNA == correctCDNA && outputProt == correctProt) ||
			strings.HasPrefix(correctCDNA, "c.-") || strings.HasPrefix(correctCDNA, "c.*") {
			continue
		}

		errorCount++
		fmt.Printf("\nWARNING: ANNOTATION MISMATCH\n")
		fmt.Printf("EXPECTED: %s|%s\n", correctCDNA, correctProt)
		fmt.Printf("RECEIVED: %s\n", annotation)

	}
	if errorCount > 7 { // from known issue
		t.Errorf("ERROR: %d variants were misannotated", errorCount)
	}
}
