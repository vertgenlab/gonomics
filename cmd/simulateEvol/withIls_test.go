package main

import (
	"fmt"
	"testing"

	"github.com/vertgenlab/gonomics/fileio"
)

var WithIlsTests = []struct {
	RootsFile              string
	TransitionMatrixFile   string
	ChromName              string
	OutPathPrefix          string
	UnitBranchLength       float64
	AncSeqFile             string
	LenSeq                 int64
	SetSeed                int64
	AncName                string
	LeafFastasOnly         bool
	SubstitutionMatrixFile string
	ExpectedPrefix         string
	ExpectedNumTopos       int64
}{
	{RootsFile: "testdata/ilsSimulate_roots.txt",
		TransitionMatrixFile:   "testdata/ilsSimulate_transMat.tsv",
		ChromName:              "test1",
		OutPathPrefix:          "testdata/ilsSimulate_out_1",
		UnitBranchLength:       0.01,
		AncSeqFile:             "",
		LenSeq:                 14,
		SetSeed:                3,
		AncName:                "Anc1",
		LeafFastasOnly:         true,
		SubstitutionMatrixFile: "",
		ExpectedPrefix:         "testdata/ilsSimulate_expected_1",
		ExpectedNumTopos:       4,
	},
	{RootsFile: "testdata/ilsSimulate_roots.txt",
		TransitionMatrixFile:   "testdata/ilsSimulate_transMat.tsv",
		ChromName:              "test2",
		OutPathPrefix:          "testdata/ilsSimulate_out_2",
		UnitBranchLength:       0.01,
		AncSeqFile:             "testdata/ilsSimulate_in_2_anc.fasta",
		LenSeq:                 50,
		SetSeed:                5,
		AncName:                "Anc2",
		LeafFastasOnly:         true,
		SubstitutionMatrixFile: "",
		ExpectedPrefix:         "testdata/ilsSimulate_expected_2",
		ExpectedNumTopos:       4,
	},
	{RootsFile: "testdata/ilsSimulate_roots.txt",
		TransitionMatrixFile:   "testdata/ilsSimulate_transMat.tsv",
		ChromName:              "test3",
		OutPathPrefix:          "testdata/ilsSimulate_out_3",
		UnitBranchLength:       .01,
		AncSeqFile:             "",
		LenSeq:                 1000,
		SetSeed:                11,
		AncName:                "Anc3",
		LeafFastasOnly:         false,
		SubstitutionMatrixFile: "",
		ExpectedPrefix:         "testdata/ilsSimulate_expected_3",
		ExpectedNumTopos:       4,
	},
}

func TestSimulateIls(t *testing.T) {
	var s IlsSettings
	for vIdx, v := range WithIlsTests {
		s = IlsSettings{
			RootsFile:              v.RootsFile,
			TransitionMatrixFile:   v.TransitionMatrixFile,
			ChromName:              v.ChromName,
			OutPathPrefix:          v.OutPathPrefix,
			UnitBranchLength:       v.UnitBranchLength,
			AncSeqFile:             v.AncSeqFile,
			LenSeq:                 v.LenSeq,
			SetSeed:                v.SetSeed,
			AncName:                v.AncName,
			LeafFastasOnly:         v.LeafFastasOnly,
			SubstitutionMatrixFile: v.SubstitutionMatrixFile,
		}

		Ils(s)

		if !fileio.AreEqual(v.OutPathPrefix+"_ils.fasta", v.ExpectedPrefix+"_ils.fasta") {
			fmt.Println(v.OutPathPrefix + "_ils.fasta")
			t.Errorf("Error in SimulateEvol ils. Output fasta %d was not as expected.", vIdx)
		}

		if !fileio.AreEqual(v.OutPathPrefix+".bed", v.ExpectedPrefix+".bed") {
			fmt.Println(v.OutPathPrefix + ".bed")
			t.Errorf("Error in SimulateEvol ils. Output bed %d was not as expected.", vIdx)
		}

		for idx := range v.ExpectedNumTopos {
			fileio.EasyRemove(fmt.Sprintf("%s_forward_evolved_topo_v%d.fasta", v.OutPathPrefix, idx))
		}

		fileio.EasyRemove(v.OutPathPrefix + ".bed")
		fileio.EasyRemove(v.OutPathPrefix + "_ils.fasta")
		fileio.EasyRemove(v.OutPathPrefix + "_anc.fasta")

	}
}
