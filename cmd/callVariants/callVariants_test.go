package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"os"
	"testing"
)

// TODO add better testing once sam simulation code is done
func TestCallVariants(t *testing.T) {
	var err error
	ref := "testdata/human_chrM.fasta"
	exp := []string{"testdata/human_chrM.bam"}
	norm := []string{"testdata/human_chrM2.bam"}
	outFile := "testdata/test_output.vcf"
	expectedFile := "testdata/test_expected.vcf"
	maxP := 1.1
	minAf := 0.0
	maxAf := 1.0
	maxStrandBias := 1.0
	minCoverage := 0
	minMapQ := 0
	minAltReads := 0
	callVariants(exp, norm, ref, outFile, maxP, minAf, maxAf, maxStrandBias, minCoverage, minMapQ, minAltReads, 1)

	if !fileio.AreEqualIgnoreComments(outFile, expectedFile) {
		t.Error("problem with variant calling")
	}

	if !t.Failed() {
		err = os.Remove(outFile)
		exception.PanicOnErr(err)
	}
}

// Start of testdata/human_chrM.fasta for reference
//
// IDX:			1234567890
// REF:			GATCACAGGT
// V1 (SNV):   	   G
// V2 (DEL):       ----
// V3 (INS):       |ATG
// V4 (DEL):       -------
// V5 (INS):       |C
// V6 (SNV):       T

func TestMultiAlleleComplexity(t *testing.T) {

	v := vcf.Vcf{
		Chr: "chrM",
		Pos: 4,
		Ref: "C",
		Alt: []string{
			"G",   // V1
			"4",   // V2
			"ATG", // V3
			"7",   // V4
			"C",   // V5
			"T",   // V6
		},
	}

	expected := vcf.Vcf{
		Chr: "chrM",
		Pos: 3,
		Ref: "TCACAGGT",
		Alt: []string{
			"TGACAGGT",
			"TGGT",
			"TCATGACAGGT",
			"T",
			"TCCACAGGT",
			"TTACAGGT",
		},
	}

	altTypes := []variantType{singleNucleotide, deletion, insertion, deletion, insertion, singleNucleotide}
	deletionIndexes := []int{1, 3}
	ref := fasta.NewSeeker("testdata/human_chrM.fasta", "")

	v = adjustAlts(v, deletionIndexes, altTypes, ref)

	if v.Chr != expected.Chr || v.Pos != expected.Pos || v.Ref != expected.Ref || len(v.Alt) != len(expected.Alt) {
		fmt.Printf("Actual:\t\t%sExpected:\t%s\n", v, expected)
		t.Error("problem with multi-allelic vcf writing")
	}

	for i := range v.Alt {
		if v.Alt[i] != expected.Alt[i] {
			fmt.Printf("Actual:\t\t%sExpected:\t%s\n", v, expected)
			t.Error("problem with multi-allelic vcf writing")
		}
	}

	err := ref.Close()
	exception.PanicOnErr(err)
}
