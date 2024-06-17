package pFasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var FaToPfaTests = []struct {
	Input          []fasta.Fasta
	Start          int
	End            int
	Chrom          string
	ExpectedPfasta PFasta
}{{Input: []fasta.Fasta{fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("TCATGACCAG")}},
	Start: 0,
	End:   9,
	Chrom: "chr1",
	ExpectedPfasta: PFasta{Name: "chr1",
		Seq: []pDna.Float32Base{
			{A: 0, C: 0, G: 0, T: 1},
			{A: 0, C: 1, G: 0, T: 0},
			{A: 1, C: 0, G: 0, T: 0},
			{A: 0, C: 0, G: 0, T: 1},
			{A: 0, C: 0, G: 1, T: 0},
			{A: 0, C: 0, G: 0, T: 0},
			{A: 1, C: 0, G: 0, T: 0},
			{A: 0, C: 1, G: 0, T: 0},
			{A: 0, C: 1, G: 0, T: 0},
			{A: 1, C: 0, G: 0, T: 0},
			{A: 0, C: 0, G: 1, T: 0},
		}}},
	{Input: []fasta.Fasta{fasta.Fasta{Name: "chr2", Seq: dna.StringToBases("TACTCCCG")},
		fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("CGACGACTTTATCC")}},
		Start: 0,
		End:   6,
		Chrom: "chr1",
		ExpectedPfasta: PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				{A: 0, C: 1, G: 0, T: 0},
				{A: 0, C: 0, G: 1, T: 0},
				{A: 1, C: 0, G: 0, T: 0},
				{A: 0, C: 1, G: 0, T: 0},
				{A: 0, C: 0, G: 1, T: 0},
				{A: 1, C: 0, G: 0, T: 0},
				{A: 0, C: 0, G: 0, T: 1},
				{A: 0, C: 0, G: 0, T: 1},
			}}},
	{Input: []fasta.Fasta{fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("TACTCCCG")},
		fasta.Fasta{Name: "chr2", Seq: dna.StringToBases("CGACGACTTTATCC")}},
		Start: 3,
		End:   7,
		Chrom: "chr1",
		ExpectedPfasta: PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				{A: 0, C: 1, G: 1, T: 1},
				{A: 0, C: 1, G: 0, T: 0},
				{A: 0, C: 1, G: 0, T: 1},
				{A: 0, C: 1, G: 0, T: 1},
				{A: 0, C: 0, G: 1, T: 1},
			}}},
}

func TestFaToPfa(t *testing.T) {
	for idx, testCase := range FaToPfaTests {
		filename := fmt.Sprintf("testdata_tools/test_faToPfa_input_%d.fa", idx)
		fasta.Write(filename, testCase.Input)

		testOutput = MultiFaToPfa()

		// faToPfa then sample(pfa) should return input fa
	}
}
