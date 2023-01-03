package motif

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var ScoreWindowTests = []struct {
	MatrixFile     string
	Seq            [][]dna.Base
	AlnStart       []int
	ExpectedScores [][]float64
	ExpectedBools  [][]bool
}{
	{MatrixFile: "testdata/jaspar.small.txt",
		Seq:            [][]dna.Base{dna.StringToBases("GCGCAGGGCAGGGCGCAGTTCAGG"), dna.StringToBases("ATGAGTTCAAGGTCAGCATGAGTTCATTGTCAGC")},
		AlnStart:       []int{0, 4, 12, 30},
		ExpectedScores: [][]float64{[]float64{43804, 21001, 34109, -1}, []float64{54010, 13722, 16827, -1}},
		ExpectedBools:  [][]bool{[]bool{true, true, true, false}, []bool{true, true, true, false}}},
}

func TestScoreWindow(t *testing.T) {
	var i, j int
	var motifs []PositionMatrix
	var currScore float64
	var currBool bool
	for _, v := range ScoreWindowTests {
		motifs = Read(v.MatrixFile, "Frequency")
		for i = range motifs {
			for j = range v.AlnStart {
				currScore, currBool = ScoreWindow(motifs[i], v.Seq[i], v.AlnStart[j])
				if currScore != v.ExpectedScores[i][j] {
					t.Errorf("Error in ScoreWindow. Score was not as expected. Expected: %v. Observed: %v.", v.ExpectedScores[i][j], currScore)
				}
				if currBool != v.ExpectedBools[i][j] {
					t.Errorf("Error in ScoreWindow. Bool was not as expected. Expected: %v. Observed: %v.", v.ExpectedBools[i][j], currBool)
				}
			}
		}
	}
}

var MatchCompTests = []struct {
	MotifFile          string
	FastaFile          string
	ChromName          string
	PropMatch          float64
	OutFile            string
	RefStart           int
	OutputAsProportion bool
	ExpectedFile       string
}{
	{MotifFile: "testdata/jaspar.vertebrate.txt",
		FastaFile: "testdata/STR012.fa",
	ChromName: "chr9",
	PropMatch: 0.95,
	OutFile: "testdata/tmp.MatchComp.bed",
	RefStart: 113944,
	OutputAsProportion: false,
	ExpectedFile: "testdata/expected.MatchComp.bed"},
	{MotifFile: "testdata/jaspar.vertebrate.txt",
		FastaFile: "testdata/STR012.fa",
		ChromName: "chr9",
		PropMatch: 0.95,
		OutFile: "testdata/tmp.OutputAsProp.MatchComp.bed",
		RefStart: 113944,
		OutputAsProportion: true,
		ExpectedFile: "testdata/expected.OutputAsProp.MatchComp.bed"},
}

func TestMatchComp(t *testing.T) {
	var err error
	var motifs []PositionMatrix
	var records []fasta.Fasta
	for _, v := range MatchCompTests {
		motifs = Read(v.MotifFile, "Frequency")
		records = fasta.Read(v.FastaFile)
		fasta.AllToUpper(records)
		MatchComp(motifs, records, v.ChromName, v.PropMatch, v.OutFile, v.RefStart, v.OutputAsProportion)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in MatchComp. Output is not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
