package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var globalAlignmentAnchorTests = []struct {
	in_maf                                   string
	species1                                 string
	species2                                 string
	species1_genome                          string
	species2_genome                          string
	gapSizeProductLimit                      int
	chrMap_filename                          string
	out_filename_prefix                      string
	diagonal                                 bool
	out_maf_expected                         string
	species1_match_bed_expected              string
	species2_match_bed_expected              string
	species1_gap_bed_expected                string
	species2_gap_bed_expected                string
	species1_gap_doNotCalculate_bed_expected string
	species2_gap_doNotCalculate_bed_expected string
	out_alignment_expected                   string
	out_species1_alignment_expected          string
	out_species2_alignment_expected          string
}{
	{"testdata/in_hg38_vs_rheMac10.toy.maf", "hg38", "rheMac10", "testdata/hg38.toy.fa", "testdata/rheMac10.toy.fa", 10000000000, "testdata/hg38_vs_rheMac10_chrMap.txt", "testdata/out_1", true, "testdata/in_hg38_vs_rheMac10.toy.filtered.12.expected.maf", "testdata/out_hg38_match.12.expected.bed", "testdata/out_rheMac10_match.12.expected.bed", "testdata/out_hg38_gap.1.expected.bed", "testdata/out_rheMac10_gap.1.expected.bed", "testdata/out_hg38_gap_doNotCalculate.1.expected.bed", "testdata/out_rheMac10_gap_doNotCalculate.1.expected.bed", "testdata/out_alignment.1.expected.tsv", "testdata/out_hg38_alignment.1.expected.bed", "testdata/out_rheMac10_alignment.1.expected.bed"},
	{"testdata/in_hg38_vs_rheMac10.toy.maf", "hg38", "rheMac10", "testdata/hg38.toy.fa", "testdata/rheMac10.toy.fa", 99, "testdata/hg38_vs_rheMac10_chrMap.txt", "testdata/out_2", true, "testdata/in_hg38_vs_rheMac10.toy.filtered.12.expected.maf", "testdata/out_hg38_match.12.expected.bed", "testdata/out_rheMac10_match.12.expected.bed", "testdata/out_hg38_gap.2.expected.bed", "testdata/out_rheMac10_gap.2.expected.bed", "testdata/out_hg38_gap_doNotCalculate.2.expected.bed", "testdata/out_rheMac10_gap_doNotCalculate.2.expected.bed", "testdata/out_alignment.2.expected.tsv", "testdata/out_hg38_alignment.2.expected.bed", "testdata/out_rheMac10_alignment.2.expected.bed"},
	{"", "hg38", "rheMac10", "testdata/hg38.toy.fa", "testdata/rheMac10.toy.fa", 10000000000, "", "testdata/out_3", true, "", "testdata/out_hg38_match.3.expected.bed", "testdata/out_rheMac10_match.3.expected.bed", "testdata/out_hg38_gap.3.expected.bed", "testdata/out_rheMac10_gap.3.expected.bed", "testdata/out_hg38_gap_doNotCalculate.3.expected.bed", "testdata/out_rheMac10_gap_doNotCalculate.3.expected.bed", "testdata/out_alignment.3.expected.tsv", "testdata/out_hg38_alignment.3.expected.bed", "testdata/out_rheMac10_alignment.3.expected.bed"}, // this test is not on the entire globalAlignmentAnchor pipeline, but only tests the helper functions matchToGap and gapToAlignment
}

func TestGlobalAlignmentAnchorTests(t *testing.T) {
	var out_maf, species1_match, species2_match, species1_gap, species2_gap, species1_gap_doNotCalculate, species2_gap_doNotCalculate, out_alignment, out_species1_alignment, out_species2_alignment string
	var err error

	for i, test := range globalAlignmentAnchorTests {

		if test.in_maf == "" {

			species1_gap_doNotCalculate = test.out_filename_prefix + "_" + test.species1 + "_gap_doNotCalculate.bed"
			species2_gap_doNotCalculate = test.out_filename_prefix + "_" + test.species2 + "_gap_doNotCalculate.bed"
			out_alignment = test.out_filename_prefix + ".alignment.tsv"
			out_species1_alignment = test.out_filename_prefix + "_" + test.species1 + "_alignment.bed"
			out_species2_alignment = test.out_filename_prefix + "_" + test.species2 + "_alignment.bed"

			species1_gap, species2_gap = matchToGap(test.species1_match_bed_expected, test.species2_match_bed_expected, test.species1_genome, test.species2_genome, test.species1, test.species2, test.gapSizeProductLimit, test.out_filename_prefix)
			gapToAlignment(species1_gap, species2_gap, test.species1_genome, test.species2_genome, test.species1, test.species2, test.out_filename_prefix)

		} else {

			out_maf = test.out_filename_prefix + ".filtered.maf"
			species1_match = test.out_filename_prefix + "_" + test.species1 + "_match.bed"
			species2_match = test.out_filename_prefix + "_" + test.species2 + "_match.bed"
			species1_gap = test.out_filename_prefix + "_" + test.species1 + "_gap.bed"
			species2_gap = test.out_filename_prefix + "_" + test.species2 + "_gap.bed"
			species1_gap_doNotCalculate = test.out_filename_prefix + "_" + test.species1 + "_gap_doNotCalculate.bed"
			species2_gap_doNotCalculate = test.out_filename_prefix + "_" + test.species2 + "_gap_doNotCalculate.bed"
			out_alignment = test.out_filename_prefix + ".alignment.tsv"
			out_species1_alignment = test.out_filename_prefix + "_" + test.species1 + "_alignment.bed"
			out_species2_alignment = test.out_filename_prefix + "_" + test.species2 + "_alignment.bed"

			globalAlignmentAnchor(test.in_maf, test.species1, test.species2, test.species1_genome, test.species2_genome, test.gapSizeProductLimit, test.chrMap_filename, test.out_filename_prefix, test.diagonal)

			if !fileio.AreEqual(out_maf, test.out_maf_expected) {
				t.Errorf("Error in out_maf, test case index: %v\n", i)
			} else {
				err = os.Remove(out_maf)
			}

			if !fileio.AreEqual(species1_match, test.species1_match_bed_expected) {
				t.Errorf("Error in species1_match_bed, test case index: %v\n", i)
			} else {
				err = os.Remove(species1_match)
			}

			if !fileio.AreEqual(species2_match, test.species2_match_bed_expected) {
				t.Errorf("Error in species2_match_bed, test case index: %v\n", i)
			} else {
				err = os.Remove(species2_match)
			}
		}

		if !fileio.AreEqual(species1_gap, test.species1_gap_bed_expected) {
			t.Errorf("Error in species1_gap_bed, test case index: %v\n", i)
		} else {
			err = os.Remove(species1_gap)
		}

		if !fileio.AreEqual(species2_gap, test.species2_gap_bed_expected) {
			t.Errorf("Error in species2_gap_bed, test case index: %v\n", i)
		} else {
			err = os.Remove(species2_gap)
		}

		if !fileio.AreEqual(species1_gap_doNotCalculate, test.species1_gap_doNotCalculate_bed_expected) {
			t.Errorf("Error in species1_gap_doNotCalculate_bed, test case index: %v\n", i)
		} else {
			err = os.Remove(species1_gap_doNotCalculate)
		}

		if !fileio.AreEqual(species2_gap_doNotCalculate, test.species2_gap_doNotCalculate_bed_expected) {
			t.Errorf("Error in species2_gap_doNotCalculate_bed, test case index: %v\n", i)
		} else {
			err = os.Remove(species2_gap_doNotCalculate)
		}

		if !fileio.AreEqual(out_alignment, test.out_alignment_expected) {
			t.Errorf("Error in out_alignment, test case index: %v\n", i)
		} else {
			err = os.Remove(out_alignment)
		}

		if !fileio.AreEqual(out_species1_alignment, test.out_species1_alignment_expected) {
			t.Errorf("Error in out_species1_alignment, test case index: %v\n", i)
		} else {
			err = os.Remove(out_species1_alignment)
		}

		if !fileio.AreEqual(out_species2_alignment, test.out_species2_alignment_expected) {
			t.Errorf("Error in out_species2_alignment, test case index: %v\n", i)
		} else {
			err = os.Remove(out_species2_alignment)
		}

		exception.PanicOnErr(err)
	}
}
