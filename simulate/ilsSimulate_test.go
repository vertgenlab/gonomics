package simulate

import (
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

func TestIlsSimulate(t *testing.T) {
	var prefix string
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
	}
}
