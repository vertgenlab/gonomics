package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var faToPfaTests = []struct {
	InFile       string
	OutDir        string
	Start      int
	End        int
	Chrom          string
}{
	{
		InFile: "testdata/test_faToPfa_input_0.fa",
		OutDir: "testdata/test_faToPfa_output_0.pfa",
		Start: 0,
		End: 10,
		Chrom: "chr1",
	}, {
		InFile: "testdata/test_faToPfa_input_0.fa",
		OutDir: "testdata/test_faToPfa_output_0_no_chr.pfa",
		Start: 0,
		End: 10,
		Chrom: "",
	}, {
		InFile: "testdata/test_faToPfa_input_1.fa",
		OutDir: "testdata/test_faToPfa_output_1.pfa",
		Start: 0,
		End: -1,
		Chrom: "chr1",
	}, {
		InFile: "testdata/test_faToPfa_input_1.fa",
		OutDir: "testdata/test_faToPfa_output_1.pfa",
		Start: 0,
		End: -1,
		Chrom: "chr1",
	},
}

func TestFaToPfa(t *testing.T) {
	var err error
	var s FaToPfaSettings
	for idx, v := range faToPfaTests {
		s = FaToPfaSettings{
			InFile: v.InFile,
			OutDir: v.OutDir,
			Start: v.Start,
			End: v.End,
			Chrom: v.Chrom,
		}

		pFaFaToPfa(s)

		testInput := fasta.Read(v.InFile)

		sampleChrom := v.Chrom
		if sampleChrom == "" && len(testInput) == 1 {
			sampleChrom = testInput[0].Name
		}


		testSample := pFasta.Sample(pFasta.Read(v.OutDir), sampleChrom)
		sampleOutputFilename := fmt.Sprintf("testdata/output_fa_%v.fa", idx)
		fasta.Write(sampleOutputFilename, []fasta.Fasta{testSample})

		testTrue := false
		for _, seq := range testInput {
			end := v.End
			if end == -1 {
				end = len(seq.Seq)
			}
			extractedSeq := fasta.Extract(seq, v.Start, end, seq.Name)
			if fasta.IsEqual(testSample, extractedSeq) {
				testTrue = true
				err = os.Remove(v.OutDir)
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