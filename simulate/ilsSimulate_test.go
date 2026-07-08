package simulate

import (
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
>>>>>>> b1a8048e (make Simulate a wrapper for funtion that only deals with data structure)
	"encoding/csv"
=======
>>>>>>> f7180d0e (test cases for ilsSimulate and matrix)
	"fmt"
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
	GenePred       string
	Seed           int64
	ExpectedPrefix string
	Precision      float64
}{
	{TransMat: "testdata/ilsSimulate_transMat.tsv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         14,
		OutName:        "test1",
		GenePred:       "testdata/debug.gp",
		Seed:           3,
		ExpectedPrefix: "testdata/ilsSimulate_expected_1",
		Precision:      1e-3,
	},
	{TransMat: "testdata/ilsSimulate_transMat.tsv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         50,
		OutName:        "test2",
		GenePred:       "testdata/debug.gp",
		Seed:           5,
		ExpectedPrefix: "testdata/ilsSimulate_expected_2",
		Precision:      1e-3,
	},
}

<<<<<<< HEAD
<<<<<<< HEAD
/// if helper functions called by main function that is used by testing, then covered on helper function

// do I need tests for the helper functions?

>>>>>>> 661b080c (skeleton code for ils simulation)
=======
>>>>>>> d4a18aea (testcase for ilsSimulate)
// probably should move this somewhere as a helper function, probably when command written
// ReadMatrix reads a tab-delimited file
// this might exist already in gonum, if not, put in numbers/matrix
// as opposed to sparse -- indices + values of non-0, implicit 0s
<<<<<<< HEAD
<<<<<<< HEAD
func readDenseFromCSV(filePath string) (*mat.Dense, error) {
=======
func readDenseFromCSV(filePath string) ([][]float64, error) {
>>>>>>> 661b080c (skeleton code for ils simulation)
=======
func readDenseFromCSV(filePath string) (*mat.Dense, error) {
>>>>>>> d4a18aea (testcase for ilsSimulate)
	// rename variables and delete this comment
	// Source - https://stackoverflow.com/a/58841827
	// Posted by SyntaxRules
	// Retrieved 2026-04-16, License - CC BY-SA 4.0
<<<<<<< HEAD
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
=======
	f, err := os.Open(filePath)
>>>>>>> d4a18aea (testcase for ilsSimulate)
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

	for i, row := range records {
		if len(row) != cols {
			return nil, fmt.Errorf("row %d has %d columns, expected %d", i, len(row), cols)
		}

		for j, cell := range row {
			val, err := strconv.ParseFloat(cell, 64)
			if err != nil {
				return nil, fmt.Errorf("parse float at row %d, col %d: %w", i, j, err)
			}
			data = append(data, val)
		}
	}

	return mat.NewDense(len(records), cols, data), nil
}

=======
>>>>>>> f7180d0e (test cases for ilsSimulate and matrix)
func TestIlsSimulate(t *testing.T) {
	var prefix string
	var expectedIls []fasta.Fasta
	var expectedBed []bed.Bed
	for _, v := range IlsSimulateTests {
		m, err := matrix.ReadDense(v.TransMat, '\t')
		if err != nil {
			t.Fatalf("error reading transition matrix: %v", err)
		}

<<<<<<< HEAD
		transMat, err := readMatrix(v.TransMat)
		// read the roots???
		observed := IlsSimulate(v.Roots, transMat, v.Length, v.Seed, v.OutName)
		// check observed
>>>>>>> 661b080c (skeleton code for ils simulation)
=======
		roots := make([]*expandedTree.ETree, len(v.Roots))
		for i, filename := range v.Roots {

			root, err := expandedTree.ReadNewick(filename)
			if err != nil {
				t.Fatalf("error reading newick %s: %v", filename, err)
			}
			roots[i] = root
		}

		expectedIls = fasta.Read(v.ExpectedPrefix + "_ils.fasta")
		expectedBed = bed.Read(v.ExpectedPrefix + ".bed")

		anc, evolved, topoRecord, ilsEvolved := SimulateIls(roots, m, int(v.Length), v.Seed, v.OutName, 0.42, v.GenePred, false, true)

		if !fasta.AllAreEqual(ilsEvolved, expectedIls) || !bed.AllAreEqual(topoRecord, expectedBed) {
			t.Errorf("observed fasta does not match expected")
			// write if wrong
			fasta.Write(v.ExpectedPrefix+"_anc.fasta", anc)
			for idx, rec := range evolved {
				prefix = fmt.Sprintf("%s_forward_evolved_topo_v%d", v.ExpectedPrefix, idx)
				fasta.Write(prefix+".fasta", rec)
			}

			bed.Write(v.ExpectedPrefix+".bed", topoRecord)
			fasta.Write(v.ExpectedPrefix+"_ils.fasta", ilsEvolved)
		} else {

		}
	}
}

var CombineIlsSeqsTests = []struct {
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
>>>>>>> d4a18aea (testcase for ilsSimulate)
	}
}
