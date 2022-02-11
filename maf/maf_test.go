package maf

import (
	//"os"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var readWriteTests = []struct {
	filename string
}{
	{"testdata/chr22.test.maf"},
}

var srcToAssemblyAndChromTests = []struct {
	src string
	assembly_expected string
	chrom_expected string
}{
	{"hg38.chr22", "hg38", "chr22"},
}

var parseTests = []struct {
	line string
	expected Maf
}{
	{"a score=407709.000000", Maf{Score: 407709.000000, Species: nil}},
}

func TestReadWrite(t *testing.T) {
	var actual []*Maf
	for _, test := range readWriteTests {
		tempFile := test.filename + ".tmp"
		actual = Read(test.filename)
		Write(tempFile, actual)
		if !fileio.AreEqualIgnoreComments(test.filename, tempFile) {
			t.Errorf("File has changed after reading and writing.")
		}
		fileio.EasyRemove(tempFile)
	}
}

func TestSrcToAssemblyAndChrom(t *testing.T) {
	var assembly_actual string
	var chrom_actual string
	for _, test := range srcToAssemblyAndChromTests {
		assembly_actual, chrom_actual = SrcToAssemblyAndChrom(test.src)
		if assembly_actual != test.assembly_expected || chrom_actual != test.chrom_expected {
			t.Errorf("SrcToAssemblyAndChrom failed.")
		}
	}
}
