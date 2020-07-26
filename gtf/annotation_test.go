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
	f := fasta.Read("testdata/test.fa")
	newSeq := make([]dna.Base, 92198969)
	f[0].Seq = append(newSeq, f[0].Seq...)
	fasta.AllToUpper(f)
	testFasta := fasta.FastaMap(f)
	tree := GenesToIntervalTree(testGtf)
	var variant *Variant
	var words, newWords, newerWords []string
	var correctCDNA, correctProt, outputCDNA, outputProt, annotation string
	var errorCount int
	for _, val := range testVcf {
		variant, _ = VcfToVariant(val, tree, testFasta)
		words = strings.Split(val.Info, "|")
		correctCDNA = words[0]
		correctProt = words[1]
		annotation = VariantToAnnotation(variant, testFasta)
		newWords = strings.Split(annotation, "|")
		newerWords = strings.Split(newWords[3], ":")
		outputCDNA = newerWords[1]
		outputProt = newWords[4]
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