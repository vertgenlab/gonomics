package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
	"testing"
)

func TestVariantToAnnotationLarge(t *testing.T) {
	testGtf := Read("testdata/TTN.gtf")
	testVcf := vcf.Read("testdata/TTN.vcf")
	f := fasta.Read("testdata/hg38.fa")
	fasta.AllToUpper(f)
	testFasta := fasta.FastaMap(f)
	tree := GenesToIntervalTree(testGtf)
	//var err error
	var variant *Variant
	var words, newWords, newerWords []string
	var correctCDNA, correctProt, outputCDNA, outputProt, annotation string
	var errorCount int
	for _, val := range testVcf {
		//fmt.Println("NEW VARIANT")
		variant, _ = VcfToVariant(val, tree, testFasta)
		//if err != nil {
		//	log.Println(err)
		//}
		words = strings.Split(val.Info,"|")
		correctCDNA = words[0]
		correctProt = words[1]
		annotation = VariantToAnnotation(variant, testFasta)
		newWords = strings.Split(annotation, "|")
		newerWords = strings.Split(newWords[3], ":")
		outputCDNA = newerWords[1]
		outputProt = newWords[4]
		//fmt.Sprintln(newerWords)
		if (outputCDNA == correctCDNA && outputProt == correctProt) ||
			strings.HasPrefix(correctCDNA, "c.-") || strings.HasPrefix(correctCDNA, "c.*") {
			continue
		}

			errorCount++
			fmt.Printf("\nWARNING: ANNOTATION MISMATCH\n")
			fmt.Printf("EXPECTED: %s|%s\n", correctCDNA, correctProt)
			fmt.Printf("RECEIVED: %s\n", annotation)

	}
	if errorCount != 0 {
		t.Errorf("ERROR: %d variants were misannotated", errorCount)
	}
}
