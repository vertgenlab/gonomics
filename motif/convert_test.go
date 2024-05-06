package motif

import (
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var PfmToPpmTests = []struct {
	PfmFile      string
	OutputFile   string
	ExpectedFile string
	Pseudocount  float64
}{
	{"testdata/expected.jaspar.txt",
		"testdata/tmp.Ppm.txt",
		"testdata/expected.Ppm.txt",
		0.1},
}

func TestPfmSliceToPpmSlice(t *testing.T) {
	var err error
	var records []PositionMatrix
	var answer []PositionMatrix
	for _, v := range PfmToPpmTests {
		records = ReadJaspar(v.PfmFile, "Frequency")
		answer = PfmSliceToPpmSlice(records, v.Pseudocount)
		WriteJaspar(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in PfmSliceToPpmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}

var PpmTests = []struct {
	PpmFile      string
	OutputFile   string
	ExpectedFile string
	GcContent    float64
}{
	{PpmFile: "testdata/expected.Ppm.txt",
		OutputFile:   "testdata/tmp.Pwm.txt",
		ExpectedFile: "testdata/expected.Pwm.txt",
		GcContent:    0.5},
	{PpmFile: "testdata/expected.Ppm.txt",
		OutputFile:   "testdata/humGc.Pwm.txt",
		ExpectedFile: "testdata/expected.humGc.Pwm.txt",
		GcContent:    0.41,
	},
}

func TestPpmSliceToPwmSlice(t *testing.T) {
	var err error
	var records []PositionMatrix
	var answer []PositionMatrix
	for _, v := range PpmTests {
		records = ReadJaspar(v.PpmFile, "Probability")
		answer = PpmSliceToPwmSlice(records, v.GcContent)
		WriteJaspar(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in PpmSliceToPwmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}

var PwmToPpmTests = []struct {
	PwmFile      string
	OutputFile   string
	ExpectedFile string
}{
	{PwmFile: "testdata/expected.Pwm.txt",
		OutputFile:   "testdata/tmp.PwmToPpm.txt",
		ExpectedFile: "testdata/expected.PwmToPpm.txt",
	},
}

func TestPwmSliceToPpmSlice(t *testing.T) {
	var records []PositionMatrix
	var answer []PositionMatrix
	var err error
	for _, v := range PwmToPpmTests {
		records = ReadJaspar(v.PwmFile, "Weight")
		answer = PwmSliceToPpmSlice(records)
		WriteJaspar(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("ERror in PwmSliceToPpmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}

var ConsensusSequencesTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	Type         string
	TieBreak     bool
}{
	{"testdata/expected.jaspar.txt",
		"testdata/tmp.jasparPFM.consensus.fa",
		"testdata/expected.jasparPFM.consensus.fa",
		"Frequency", true},
}

func TestConsensusSequences(t *testing.T) {
	var err error
	var input []PositionMatrix
	var records []fasta.Fasta
	for _, v := range ConsensusSequencesTests {
		input = ReadJaspar(v.InFile, v.Type)
		records = ConsensusSequences(input, v.TieBreak)
		fasta.Write(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in Consensus Sequences. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var ReverseComplementTests = []struct {
	InFile              string
	RevCompOutFile      string
	RevCompExpectedFile string
}{
	{"testdata/jaspar.vertebrate.txt",
		"testdata/tmp.RevCompConsensus.txt",
		"testdata/expected.RevCompConsensus.txt"},
}

func TestReverseComplement(t *testing.T) {
	var err error
	var input []PositionMatrix
	var reversed []PositionMatrix
	var reverseSequences []fasta.Fasta
	for _, v := range ReverseComplementTests {
		input = ReadJaspar(v.InFile, "Frequency")
		reversed = ReverseComplementAll(input)
		reverseSequences = ConsensusSequences(reversed, false)
		fasta.Write(v.RevCompOutFile, reverseSequences)
		if !fileio.AreEqual(v.RevCompOutFile, v.RevCompExpectedFile) {
			t.Errorf("Error in ReverseComplement. Output was not as expected.")
		} else {
			err = os.Remove(v.RevCompOutFile)
			exception.PanicOnErr(err)
		}
	}
}
