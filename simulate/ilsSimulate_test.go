package simulate

import (
<<<<<<< HEAD
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

			root, err := expandedTree.ReadNewick(filename)
=======
=======
	"github.com/vertgenlab/gonomics/numbers/matrix"
	"github.com/vertgenlab/gonomics/fasta"
	"encoding/csv"
	"os"
	"strconv"
	"testing"
)

var IlsSimulateTests = []struct {
	TransMat    string
	Roots 	[]*expandedTree.ETree
	Length    int64
	OutName   string
	Seed	  int64
	Expected  fasta.Fasta
	Precision float64
}{
	{TransMat: "testdata/transMat_1.csv",
		Roots: "testdata/roots_1.????",
		Length: 10,
		OutName: "test1"
		Seed: 3,
		Expected: "testdata/ilsSimulate_expected_1.fasta"
		Precision: 1e-3,
	},
}

/// if helper functions called by main function that is used by testing, then covered on helper function

// do I need tests for the helper functions?

>>>>>>> 661b080c (skeleton code for ils simulation)
// probably should move this somewhere as a helper function, probably when command written
// ReadMatrix reads a tab-delimited file
// this might exist already in gonum, if not, put in numbers/matrix
// as opposed to sparse -- indices + values of non-0, implicit 0s
<<<<<<< HEAD
func readDenseFromCSV(filePath string) (*mat.Dense, error) {
=======
func readDenseFromCSV(filePath string) ([][]float64, error) {
>>>>>>> 661b080c (skeleton code for ils simulation)
	// rename variables and delete this comment
	// Source - https://stackoverflow.com/a/58841827
	// Posted by SyntaxRules
	// Retrieved 2026-04-16, License - CC BY-SA 4.0
<<<<<<< HEAD
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
=======
	ff, err := os.Open(filePath)
    if err != nil {
        log.Fatal("Unable to read input file " + filePath, err)
    }
    defer f.Close()

    csvReader := csv.NewReader(f)
	csvReader.Comma = '\t'
    records, err := csvReader.ReadAll()
    if err != nil {
        log.Fatal("Unable to parse file as CSV for " + filePath, err)
    }

	matrix := make([][]float64, len(records))

	for i := range records {
		matrix[i] = make([]float64, len(records[i]))
		for j := range records[i] {
			val, err := strconv.ParseFloat(records[i][j], 64)
			if err != nil {
				return nil, err
			}
			matrix[i][j] = val
		}
	}

	return matrix, nil
}

func TestIlsSimulate(t *testing.T) {
	var err error
	var observed [][]float64
	for _, v := range ilsSimulateTests {
		// func SimulateIls(roots []*expandedTree.ETree, transitionMat [][]float64, totalLength int, seed int, outName string) {

		transMat, err := readMatrix(v.TransMat)
		// read the roots??? 
		observed := IlsSimulate(v.Roots, transMat, v.Length, v.Seed, v.OutName)
		// check observed
>>>>>>> 661b080c (skeleton code for ils simulation)
	}
}
