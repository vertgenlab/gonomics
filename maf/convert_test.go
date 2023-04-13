package maf

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

var MafSLinesAreEqualTests = []struct {
	a        MafSLine
	b        MafSLine
	expected bool
}{
	{MafSLine{Src: "hg38.chr6_GL000254v2_alt", Start: 331, Size: 38, Strand: true, SrcSize: 4827813, Seq: dna.StringToBases("CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGTGGGGG")}, MafSLine{Src: "hg38.chr6_GL000254v2_alt", Start: 331, Size: 38, Strand: true, SrcSize: 4827813, Seq: dna.StringToBases("CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGTGGGGG")}, true},
	{MafSLine{Src: "hg38.chr6_GL000254v2_alt", Start: 331, Size: 38, Strand: true, SrcSize: 4827813, Seq: dna.StringToBases("CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGTGGGGG")}, MafSLine{Src: "hg38.chr6_GL000254v2_alt", Start: 331, Size: 38, Strand: true, SrcSize: 4827813, Seq: dna.StringToBases("CTGTAGTCTGTCAGATATGGGTGGAGTGGGGGAGGGGG")}, false},
}

var MafILinesAreEqualTests = []struct {
	a        MafILine
	b        MafILine
	expected bool
}{
	{MafILine{Src: "rheMac10.chr1", LeftStatus: 'N', LeftCount: 0, RightStatus: 'I', RightCount: 35}, MafILine{Src: "rheMac10.chr1", LeftStatus: 'N', LeftCount: 0, RightStatus: 'I', RightCount: 35}, true},
	{MafILine{Src: "rheMac10.chr1", LeftStatus: 'N', LeftCount: 0, RightStatus: 'I', RightCount: 35}, MafILine{Src: "rheMac10.chr1", LeftStatus: 'N', LeftCount: 0, RightStatus: 'I', RightCount: 33}, false},
}

var MafELinesAreEqualTests = []struct {
	a        MafELine
	b        MafELine
	expected bool
}{
	{MafELine{Src: "rheMac10.chr1", Start: 122101540, Size: 35, Strand: false, SrcSize: 223616942, Status: 'I'}, MafELine{Src: "rheMac10.chr1", Start: 122101540, Size: 35, Strand: false, SrcSize: 223616942, Status: 'I'}, true},
	{MafELine{Src: "rheMac10.chr1", Start: 122101540, Size: 35, Strand: false, SrcSize: 223616942, Status: 'I'}, MafELine{Src: "rheMac10.chr1", Start: 122101540, Size: 35, Strand: false, SrcSize: 223616942, Status: 'N'}, false},
}

func TestMafSLinesAreEqual(t *testing.T) {
	for _, test := range MafSLinesAreEqualTests {
		if MafSLinesAreEqual(test.a, test.b) != test.expected {
			t.Errorf("MafSLinesAreEqual failed.")
		}
	}
}

func TestMafILinesAreEqual(t *testing.T) {
	for _, test := range MafILinesAreEqualTests {
		if MafILinesAreEqual(test.a, test.b) != test.expected {
			t.Errorf("MafILinesAreEqual failed.")
		}
	}
}

func TestMafELinesAreEqual(t *testing.T) {
	for _, test := range MafELinesAreEqualTests {
		if MafELinesAreEqual(test.a, test.b) != test.expected {
			t.Errorf("MafELinesAreEqual failed.")
		}
	}
}
