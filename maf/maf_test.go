package maf

import (
	//"os"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

var readWriteTests = []struct {
	filename string
}{
	{"testdata/chr22.test.maf"},
}

var srcToAssemblyAndChromTests = []struct {
	src               string
	assembly_expected string
	chrom_expected    string
}{
	{"hg38.chr22", "hg38", "chr22"},
}

var parseMafALineTests = []struct {
	line     string
	expected Maf
}{
	{"a score=407709.000000", Maf{Score: 407709.000000, Species: nil}},
}

var parseMafSLineTests = []struct {
	line     string
	expected MafSLine
}{
	{"s hg38.chr6_GL000254v2_alt      331 38 +   4827813 CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGTGGGGG", MafSLine{Src: "hg38.chr6_GL000254v2_alt", Start: 331, Size: 38, Strand: true, SrcSize: 4827813, Seq: dna.StringToBases("CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGTGGGGG")}},
}

var parseMafIStatusTests = []struct {
	s        string
	expected rune
}{
	{"C", 'C'},
	{"I", 'I'},
	{"N", 'N'},
	{"n", 'n'},
	{"M", 'M'},
	{"T", 'T'},
}

var parseMafILineTests = []struct {
	line     string
	expected MafILine
}{
	{"i rheMac10.chr1               N 0 I 35", MafILine{Src: "rheMac10.chr1", LeftStatus: 'N', LeftCount: 0, RightStatus: 'I', RightCount: 35}},
}

var parseMafEStatusTests = []struct {
	s        string
	expected rune
}{
	{"C", 'C'},
	{"I", 'I'},
	{"M", 'M'},
	{"n", 'n'},
	{"T", 'T'},
}

var parseMafELineTests = []struct {
	line     string
	expected MafELine
}{
	{"e rheMac10.chr1               122101540 35 - 223616942 I", MafELine{Src: "rheMac10.chr1", Start: 122101540, Size: 35, Strand: false, SrcSize: 223616942, Status: 'I'}},
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

func TestParseMafALine(t *testing.T) {
	var aLine_actual *Maf
	for _, test := range parseMafALineTests {
		aLine_actual = parseMafALine(test.line)
		if !MafsAreEqual(*aLine_actual, test.expected) {
			t.Errorf("parseMafALine failed. &aLine_actual: %v, expected: %v", &aLine_actual, test.expected)
		}
	}
}

func TestParseMafSLine(t *testing.T) {
	var sLine_actual *MafSLine
	for _, test := range parseMafSLineTests {
		sLine_actual = parseMafSLine(test.line)
		if !MafSLinesAreEqual(*sLine_actual, test.expected) {
			t.Errorf("parseMafSLine failed. &sLine_actual: %v, expected: %v", &sLine_actual, test.expected)
		}
	}
}

func TestParseMafIStatus(t *testing.T) {
	var actual rune
	for _, test := range parseMafIStatusTests {
		actual = parseMafIStatus(test.s)
		if actual != test.expected {
			t.Errorf("parseMafIStatus failed. actual: %v, expected: %v", actual, test.expected)
		}
	}
}

func TestParseMafILine(t *testing.T) {
	var iLine_actual *MafILine
	for _, test := range parseMafILineTests {
		iLine_actual = parseMafILine(test.line)
		if !MafILinesAreEqual(*iLine_actual, test.expected) {
			t.Errorf("parseMafILine failed. &iLine_actual: %v, expected: %v", &iLine_actual, test.expected)
		}
	}
}

func TestParseMafEStatus(t *testing.T) {
	var actual rune
	for _, test := range parseMafEStatusTests {
		actual = parseMafEStatus(test.s)
		if actual != test.expected {
			t.Errorf("parseMafEStatus failed. actual: %v, expected: %v", actual, test.expected)
		}
	}
}

func TestParseMafELine(t *testing.T) {
	var eLine_actual *MafELine
	for _, test := range parseMafELineTests {
		eLine_actual = parseMafELine(test.line)
		if !MafELinesAreEqual(*eLine_actual, test.expected) {
			t.Errorf("parseMafELine failed. &eLine_actual: %v, expected: %v", &eLine_actual, test.expected)
		}
	}
}
