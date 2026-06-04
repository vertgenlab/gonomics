package simulate

import (
	"fmt"
	"path/filepath"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers/matrix"
)

var IlsSimulateTests = []struct {
	TransMat       string
	Roots          []string
	Length         int64
	OutName        string
	Seed           int64
	ExpectedPrefix string
	Precision      float64
}{
	{TransMat: "testdata/ilsSimulate_transMat.tsv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         14,
		OutName:        "test1",
		Seed:           3,
		ExpectedPrefix: "testdata/ilsSimulate_expected_1",
		Precision:      1e-3,
	},
	{TransMat: "testdata/ilsSimulate_transMat.tsv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         50,
		OutName:        "test2",
		Seed:           5,
		ExpectedPrefix: "testdata/ilsSimulate_expected_2",
		Precision:      1e-3,
	},
}

func TestIlsSimulate(t *testing.T) {
	var expectedIls []fasta.Fasta
	var expectedBed []bed.Bed
	for _, v := range IlsSimulateTests {
		m, err := matrix.ReadDense(v.TransMat, '\t')
		if err != nil {
			t.Fatalf("error reading transition matrix: %v", err)
		}

		roots := make([]*expandedTree.ETree, len(v.Roots))
		for i, filename := range v.Roots {

<<<<<<< HEAD
			root, err := expandedTree.ReadNewick(filename)
=======
// probably should move this somewhere as a helper function, probably when command written
// ReadMatrix reads a tab-delimited file
// this might exist already in gonum, if not, put in numbers/matrix
// as opposed to sparse -- indices + values of non-0, implicit 0s
func readDenseFromCSV(filePath string) (*mat.Dense, error) {
	// rename variables and delete this comment
	// Source - https://stackoverflow.com/a/58841827
	// Posted by SyntaxRules
	// Retrieved 2026-04-16, License - CC BY-SA 4.0
	f, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("open %q: %w", filePath, err)
	}
	defer f.Close()

	reader := csv.NewReader(f)
	reader.Comma = '\t'

	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("parse %q as CSV: %w", filePath, err)
	}

	if len(records) == 0 {
		return mat.NewDense(0, 0, nil), nil
	}

	cols := len(records[0])
	data := make([]float64, 0, len(records)*cols)

	for i := range records {
		matrix[i] = make([]float64, len(records[i]))
			val, err := strconv.ParseFloat(records[i][j], 64)
	forwardEvolved     [][]fasta.Fasta
	statePath          []int
	expectedIls        []fasta.Fasta
	expectedPathLength int
}{
	{forwardEvolved: [][]fasta.Fasta{
		{
			{Name: "A", Seq: dna.StringToBases("AAAA")},
			{Name: "B", Seq: dna.StringToBases("CCCC")},
		},
		{
			{Name: "B", Seq: dna.StringToBases("GGGG")},
			{Name: "A", Seq: dna.StringToBases("TTTT")}},
	},
		statePath: []int{0, 1, 1, 0},
		expectedIls: []fasta.Fasta{
			{Name: "A", Seq: dna.StringToBases("ATTA")},
			{Name: "B", Seq: dna.StringToBases("CGGC")},
		},
		expectedPathLength: 3,
	},
	{forwardEvolved: [][]fasta.Fasta{
		{
			{Name: "A", Seq: dna.StringToBases("AAAAAAAAAA")},
			{Name: "B", Seq: dna.StringToBases("CCCCCCCCCC")},
		},
		{
			{Name: "B", Seq: dna.StringToBases("GGGGGGGGGG")},
			{Name: "A", Seq: dna.StringToBases("TTTTTTTTTT")},
		},
		{
			{Name: "A", Seq: dna.StringToBases("AAACGTACGT")},
			{Name: "B", Seq: dna.StringToBases("TTTGCATGCA")},
		},
	},
		statePath: []int{0, 0, 0, 1, 2, 2, 1, 0, 2, 1},
		expectedIls: []fasta.Fasta{
			{Name: "A", Seq: dna.StringToBases("AAATGTTAGT")},
			{Name: "B", Seq: dna.StringToBases("CCCGCAGCCG")},
		},
		expectedPathLength: 7,
	},
}

<<<<<<< HEAD
func TestCombineIlsSeqsTests(t *testing.T) {
	var states []int

	for _, v := range CombineIlsSeqsTests {
		states = v.statePath

		stateRecord, outIls := CombineIlsSeqs(v.forwardEvolved, states, "ilsState")

		if !fasta.AllAreEqual(outIls, v.expectedIls) {
			t.Errorf("Expected and combined sequences are not equal")
		}

		if len(stateRecord) != v.expectedPathLength {
			t.Errorf("expected %d BED records, got %d", v.expectedPathLength, len(stateRecord))
		}
=======
		transMat, err := readMatrix(v.TransMat)
		// read the roots???
		observed := IlsSimulate(v.Roots, transMat, v.Length, v.Seed, v.OutName)
		// check observed
>>>>>>> b1a8048e (make Simulate a wrapper for funtion that only deals with data structure)
	}
}
