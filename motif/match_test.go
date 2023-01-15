package motif

import (
	"fmt"
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
		ExpectedScores: [][]float64{{43804, 21001, 34109, -1}, {54010, 13722, 16827, -1}},
		ExpectedBools:  [][]bool{{true, true, true, false}, {true, true, true, false}}},
}

func TestScoreWindow(t *testing.T) {
	var i, j int
	var motifs []PositionMatrix
	var currScore float64
	var currBool bool
	for _, v := range ScoreWindowTests {
		motifs = ReadJaspar(v.MatrixFile, "Frequency")
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
		FastaFile:          "testdata/STR012.fa",
		ChromName:          "chr9",
		PropMatch:          0.95,
		OutFile:            "testdata/tmp.MatchComp.bed",
		RefStart:           113944,
		OutputAsProportion: false,
		ExpectedFile:       "testdata/expected.MatchComp.bed"},
	{MotifFile: "testdata/jaspar.vertebrate.txt",
		FastaFile:          "testdata/STR012.fa",
		ChromName:          "chr9",
		PropMatch:          0.95,
		OutFile:            "testdata/tmp.OutputAsProp.MatchComp.bed",
		RefStart:           113944,
		OutputAsProportion: true,
		ExpectedFile:       "testdata/expected.OutputAsProp.MatchComp.bed"},
}

func TestMatchComp(t *testing.T) {
	var err error
	var motifs []PositionMatrix
	var records []fasta.Fasta
	for _, v := range MatchCompTests {
		motifs = ReadJaspar(v.MotifFile, "Frequency")
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

var RankMatrixTests = []struct {
	PwmFile      string
	OutFile      string
	ExpectedFile string
}{
	{PwmFile: "testdata/expected.Pwm.txt",
		OutFile:      "testdata/tmp.RankMatrix.txt",
		ExpectedFile: "testdata/expected.RankMatrix.txt",
	},
}

func TestRankTensors(t *testing.T) {
	var err error
	var motifs []PositionMatrix
	var i int
	var out *fileio.EasyWriter
	var currRankMatrix [][]RankTensorElement
	for _, v := range RankMatrixTests {
		out = fileio.EasyCreate(v.OutFile)
		motifs = ReadJaspar(v.PwmFile, "Weight")
		for i = range motifs {
			currRankMatrix = initializeRankTensor(motifs[i])
			_, err = fmt.Fprintf(out, rankTensorToString(currRankMatrix))
			exception.PanicOnErr(err)
		}
		err = out.Close()
		exception.PanicOnErr(err)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in RankMatrix test. Not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var BuildKmerHashTests = []struct {
	PwmFile             string
	ThresholdProportion float64
	OutFile             string
	ExpectedLengths     []int
}{
	{PwmFile: "testdata/pwm.small.txt",
		ThresholdProportion: 0.95,
		ExpectedLengths:     []int{6, 30}},
	{PwmFile: "testdata/pwm.small.txt",
		ThresholdProportion: 0.8,
		ExpectedLengths:     []int{104, 1705}},
	{PwmFile: "testdata/pwm.small.txt",
		ThresholdProportion: 0.5,
		ExpectedLengths:     []int{1658, 123496}},
}

func TestBuildKmerHash(t *testing.T) {
	var motifs []PositionMatrix
	var answer map[uint64]float64
	for _, v := range BuildKmerHashTests {
		motifs = ReadJaspar(v.PwmFile, "Weight")
		for i := range motifs {
			answer = buildKmerHash(motifs[i], v.ThresholdProportion)
			if len(answer) != v.ExpectedLengths[i] {
				t.Errorf("Error in buildKmerHash. Expected Hash Size: %d. Found: %d.", len(answer), v.ExpectedLengths[i])
			}
		}
	}
}
