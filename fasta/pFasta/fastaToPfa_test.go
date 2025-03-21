package pFasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var FaToPfaTests = []struct {
	InputFilename string
	Start         int
	End           int
	Chrom         string
}{
	{
		InputFilename: "testdata/test_faToPfa_input_0.fa",
		Start:         0,
		End:           10,
		Chrom:         "chr1",
	}, {
		InputFilename: "testdata/test_faToPfa_input_0.fa",
		Start:         0,
		End:           10,
		Chrom:         "",
	}, {
		InputFilename: "testdata/test_faToPfa_input_1.fa",
		Start:         0,
		End:           -1,
		Chrom:         "chr1",
	}, {
		InputFilename: "testdata/test_faToPfa_input_2.fa",
		Start:         3,
		End:           8,
		Chrom:         "chr1",
	},
}

func TestFaToPfa(t *testing.T) {
	var err error
	for idx, v := range FaToPfaTests {
		//running function
		testOutput := MultiFaToPfa(v.InputFilename, v.Start, v.End, v.Chrom)
		outputFilename := fmt.Sprintf("testdata/output_%v.pfa", idx)
		Write(outputFilename, []PFasta{testOutput})

		testInput := fasta.Read(v.InputFilename)

		// sample to get fasta from pfasta
		sampleChrom := v.Chrom
		if sampleChrom == "" && len(testInput) == 1 {
			sampleChrom = testInput[0].Name
		}

		// testOutput is 1-hot pFa, sampling from testOutput should return (subsequence of) input Fasta
		testSample := Sample([]PFasta{testOutput}, sampleChrom)
		sampleOutputFilename := fmt.Sprintf("testdata/output_fa_%v.fa", idx)
		fasta.Write(sampleOutputFilename, []fasta.Fasta{testSample})

		// compare all sequences in input fasta to testOutput
		testTrue := false

		for _, seq := range testInput {
			end := v.End
			if end == -1 {
				end = len(seq.Seq)
			}

			extractedSeq := fasta.Extract(seq, v.Start, end, seq.Name)
			if fasta.IsEqual(testSample, extractedSeq) {
				testTrue = true
				err = os.Remove(outputFilename)
				exception.PanicOnErr(err)
				err = os.Remove(sampleOutputFilename)
				exception.PanicOnErr(err)
			}
		}

		if !testTrue {
			t.Errorf("Error: in faToPfa. Output not as expected.\n")
		}
	}
}
