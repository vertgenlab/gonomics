package motif

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var PfmToPpmTests = []struct {
	PfmFile string
	OutputFile string
	ExpectedFile string
	Pseudocount float64
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
		records = Read(v.PfmFile, "Frequency")
		answer = PfmSliceToPpmSlice(records, v.Pseudocount)
		Write(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in PfmSliceToPpmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}

var PpmTests = []struct {
	PpmFile string
	OutputFile string
	ExpectedFile string
}{
	{"testdata/expected.Ppm.txt",
		"testdata/tmp.Pwm.txt",
	"testdata/expected.Pwm.txt"},
}

func TestPpmSliceToPwmSlice(t *testing.T) {
	var err error
	var records []PositionMatrix
	var answer []PositionMatrix
	for _, v := range PpmTests {
		records = Read(v.PpmFile, "Probability")
		answer = PpmSliceToPwmSlice(records)
		Write(v.OutputFile, answer)
		if !fileio.AreEqual(v.OutputFile, v.ExpectedFile) {
			t.Errorf("Error in PpmSliceToPwmSlice. Output was not as expected.")
		} else {
			err = os.Remove(v.OutputFile)
			exception.PanicOnErr(err)
		}
	}
}

var ConsensusSequencesTests = []struct {
	InFile string
	OutFile string
	ExpectedFile string
	Type string
}{
	{"testdata/expected.jaspar.txt",
		"testdata/tmp.jasparPFM.consensus.fa",
	"testdata/expected.jasparPFM.consensus.fa",
	"Frequency"},
}

func TestConsensusSequences(t *testing.T) {
	var err error
	var input []PositionMatrix
	var records []fasta.Fasta
	for _, v := range ConsensusSequencesTests {
		input = Read(v.InFile, v.Type)
		records = ConsensusSequences(input)
		fasta.Write(v.OutFile, records)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("ERror in Concensus Sequences. Output not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
