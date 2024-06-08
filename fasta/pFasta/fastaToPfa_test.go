package pFasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fasta"
)

var faToPfaTests = []struct {
	Input          fasta.Fasta
	Start          int
	End            int
	ExpectedPfasta PFasta
}{{Input: fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("TCATGACCAG")},
	Start: 0,
	End:   9,
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
	{Input: fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("CGACGACTTTATCC")},
		Start: 0,
		End:   6,
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
	{Input: fasta.Fasta{Name: "chr1", Seq: dna.StringToBases("TACTCCCG")},
		Start: 3,
		End:   7,
		ExpectedPfasta: PFasta{Name: "chr1",
			Seq: []pDna.Float32Base{
				{A: 0, C: 1, G: 1, T: 1},
				{A: 0, C: 1, G: 0, T: 0},
				{A: 0, C: 1, G: 0, T: 1},
				{A: 0, C: 1, G: 0, T: 1},
				{A: 0, C: 0, G: 1, T: 1},
			}}},
}
