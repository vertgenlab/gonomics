package simulate

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/expandedTree"
	"github.com/vertgenlab/gonomics/fasta"
	"gonum.org/v1/gonum/mat"
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
	{TransMat: "testdata/ilsSimulate_transMat.csv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         14,
		OutName:        "test1",
		GenePred:       "testdata/debug.gp",
		Seed:           3,
		ExpectedPrefix: "testdata/ilsSimulate_expected_1",
		Precision:      1e-3,
	},
	{TransMat: "testdata/ilsSimulate_transMat.csv",
		Roots:          []string{"testdata/ilsSimulate_v0.nh", "testdata/ilsSimulate_v1.nh", "testdata/ilsSimulate_v2.nh", "testdata/ilsSimulate_v3.nh"},
		Length:         50,
		OutName:        "test2",
		GenePred:       "testdata/debug.gp",
		Seed:           5,
		ExpectedPrefix: "testdata/ilsSimulate_expected_2",
		Precision:      1e-3,
	},
}

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

func TestIlsSimulate(t *testing.T) {
	var prefix string
	var expectedIls []fasta.Fasta
	var expectedBed []bed.Bed
	for _, v := range IlsSimulateTests {
		transMat, err := readDenseFromCSV(v.TransMat)
		if err != nil {
			t.Fatalf("error reading transition matrix: %v", err)
		}

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

		anc, evolved, topoRecord, ilsEvolved := SimulateIls(roots, transMat, int(v.Length), v.Seed, v.OutName, 0.42, v.GenePred, false, true)

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
			{Name: "A", Seq: dna.StringToBases("AAAAAAAA")},
			{Name: "B", Seq: dna.StringToBases("CCCCCCCC")},
		},
		{
			{Name: "B", Seq: dna.StringToBases("GGGGGGGG")},
			{Name: "A", Seq: dna.StringToBases("TTTTTTTT")},
		},
		{
			{Name: "A", Seq: dna.StringToBases("ACGTACGT")},
			{Name: "B", Seq: dna.StringToBases("TGCATGCA")},
		},
	},
		statePath: []int{0, 1, 2, 2, 1, 0, 2, 1},
		expectedIls: []fasta.Fasta{
			{Name: "A", Seq: dna.StringToBases("ATGTTAGT")},
			{Name: "B", Seq: dna.StringToBases("CGCAGCCG")},
		},
		expectedPathLength: 7,
	},
}

func TestCombineIlsSeqsTests(t *testing.T) {
	var states []int

	for _, v := range CombineIlsSeqsTests {
		states = v.statePath

		stateRecord, outIls := CombineIlsSeqs(v.forwardEvolved, states, "ilsState")
		// for _, seq := range outIls {
		// 	fmt.Println(dna.BasesToString(seq.Seq))
		// }

		// for idx := range outIls {
		// 	if dna.BasesToString(outIls[idx].Seq) != dna.BasesToString(v.expectedIls[idx].Seq) {
		// 		t.Errorf("expected %s , got %s", dna.BasesToString(v.expectedIls[idx].Seq), dna.BasesToString(outIls[idx].Seq))
		// 	}
		// }

		/// OR I could just do a generic "these seqs don't match up"
		if !fasta.AllAreEqual(outIls, v.expectedIls) {
			t.Errorf("Expected and combined sequences are not equal")
		}

		if len(stateRecord) != v.expectedPathLength {
			t.Errorf("expected %d BED records, got %d", v.expectedPathLength, len(stateRecord))
		}
	}
}
