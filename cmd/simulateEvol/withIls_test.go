package main

import (
	"fmt"
	"testing"

	"github.com/vertgenlab/gonomics/fasta"
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
	LeafFastasOnly         bool
	SubstitutionMatrixFile string
	ExpectedPrefix         string
}{
	{RootsFile: "testdata/ilsSimulate_roots.txt",
		TransitionMatrixFile:   "testdata/ilsSimulate_transMat.tsv",
		ChromName:              "test1",
		OutPathPrefix:          "testdata/ilsSimulate_out_1",
		UnitBranchLength:       1.0,
		AncSeqFile:             "",
		LenSeq:                 14,
		SetSeed:                3,
		LeafFastasOnly:         true,
		SubstitutionMatrixFile: "",
		ExpectedPrefix:         "testdata/ilsSimulate_expected_1",
	},
	{RootsFile: "testdata/ilsSimulate_roots.txt",
		TransitionMatrixFile:   "testdata/ilsSimulate_transMat.tsv",
		ChromName:              "test2",
		OutPathPrefix:          "testdata/ilsSimulate_out_2",
		UnitBranchLength:       1.0,
		AncSeqFile:             "testdata/ilsSimulate_expected_2_anc.fasta",
		LenSeq:                 50,
		SetSeed:                5,
		LeafFastasOnly:         true,
		SubstitutionMatrixFile: "",
		ExpectedPrefix:         "testdata/ilsSimulate_expected_2",
	},
}

func TestSimulateIls(t *testing.T) {
	var s IlsSettings
	var numTopos int
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

		numTopos = len(fasta.Read(v.OutPathPrefix + "_ils.fasta"))

		for idx := range numTopos {
			fileio.EasyRemove(fmt.Sprintf("%s_forward_evolved_topo_v%d.fasta", v.OutPathPrefix, idx))
		}

		fileio.EasyRemove(v.OutPathPrefix + ".bed")
		fileio.EasyRemove(v.OutPathPrefix + "_ils.fasta")
		fileio.EasyRemove(v.OutPathPrefix + "_anc.fasta")

	}
}
