package simulate

import (
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

// probably should move this somewhere as a helper function, probably when command written
// ReadMatrix reads a tab-delimited file
// this might exist already in gonum, if not, put in numbers/matrix
// as opposed to sparse -- indices + values of non-0, implicit 0s
func readDenseFromCSV(filePath string) ([][]float64, error) {
	// rename variables and delete this comment
	// Source - https://stackoverflow.com/a/58841827
	// Posted by SyntaxRules
	// Retrieved 2026-04-16, License - CC BY-SA 4.0
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
	}
}
