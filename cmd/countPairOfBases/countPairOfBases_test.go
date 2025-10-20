package main

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

var CountPairOfBasesTest = []struct {
	fastaFile   string
	chromName   string
	baseone     string
	basetwo     string
	outfile     string
	bedFile     string
	compare     bool
	expectedTxt string
}{
	{"./testdata/singlegenome/chr1.fa", "chr1", "C", "G", "single_genome_chr1out.bed", "", false, "./testdata/single_genome_chr1expected.txt"},
	{"./testdata/singlegenome/chr11.fa", "chr11", "C", "G", "single_genome_chr11out.bed", "", false, "./testdata/single_genome_chr11expected.txt"},
	{"./testdata/singlegenome/chr1.fa", "chr1", "C", "G", "single_genome_chr1bedOut.bed", "./testdata/single_genome_chr1test.bed", false, "./testdata/single_genome_chr1bedExpected.bed"},
	{"./testdata/singlegenome/chr11.fa", "chr11", "C", "G", "single_genome_chr11bedOut.bed", "./testdata/single_genome_chr11test.bed", false, "./testdata/single_genome_chr11bedExpected.bed"},
	{"./testdata/twogenome/chr1.fa", "chr1", "C", "G", "two_genome_chr1out.bed", "", true, "./testdata/two_genome_chr1expected.txt"},
	{"./testdata/twogenome/chr14.fa", "chr14", "C", "G", "two_genome_chr14out.bed", "", true, "./testdata/two_genome_chr14expected.txt"},
	{"./testdata/twogenome/chr1.fa", "chr1", "C", "G", "two_genome_chr1bedOut.bed", "./testdata/two_genome_chr1test.bed", true, "./testdata/two_genome_chr1bedExpected.txt"},
	{"./testdata/twogenome/chr14.fa", "chr14", "C", "G", "two_genome_chr14bedOut.bed", "./testdata/two_genome_chr14test.bed", true, "./testdata/two_genome_chr14bedExpected.txt"},
}

func TestCountPairOfBases(t *testing.T) {
	var err error
	var s Settings
	for _, v := range CountPairOfBasesTest {
		s = Settings{
			InFa:    v.fastaFile,
			Chrom:   v.chromName,
			BaseOne: v.baseone,
			BaseTwo: v.basetwo,
			Outfile: v.outfile,
			BedFile: v.bedFile,
			Compare: v.compare,
		}
		countPairOfBases(s)

		if !fileio.AreEqual(v.expectedTxt, v.outfile) {
			t.Errorf("Error in countPairOfBases.go, expected output file != actual output file")
		} else {
			err = os.Remove(v.outfile)
			exception.PanicOnErr(err)
		}
	}
}
