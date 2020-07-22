package gtf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
	"testing"
)

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

		fmt.Printf("\nWARNING: ANNOTATION MISMATCH\n")
		fmt.Printf("EXPECTED: %s\n", correctString)
		fmt.Printf("RECEIVED: %s\n", annotation)

	}
	if errorCount != 0 {
		t.Errorf("ERROR: %d variants were misannotated", errorCount)
	}
}
