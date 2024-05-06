package motif

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"os"
	"testing"
)

var MatchCompTests = []struct {
	MotifFile          string
	MotifType          string
	FastaFile          string
	ChromName          string
	PropMatch          float64
	OutFile            string
	RefStart           int
	ResidualWindowSize int
	OutputAsProportion bool
	Pseudocounts       float64
	ExpectedFile       string
	EnforceStrandMatch bool
	GcContent          float64
}{
	{MotifFile: "testdata/myMotifFile.txt",
		MotifType:          "Frequency",
		FastaFile:          "testdata/myAln.fa",
		ChromName:          "chr1",
		PropMatch:          0.1,
		OutFile:            "testdata/tmp.myAln.bed",
		RefStart:           100,
		OutputAsProportion: true,
		ResidualWindowSize: 1,
		Pseudocounts:       0.1,
		ExpectedFile:       "testdata/expected.myAln.bed",
		EnforceStrandMatch: false,
		GcContent:          0.5,
	},
	{MotifFile: "testdata/myMotifFile.txt",
		MotifType:          "Frequency",
		FastaFile:          "testdata/myAln.fa",
		ChromName:          "chr1",
		PropMatch:          0.1,
		OutFile:            "testdata/tmp.myAln.highPseudo.bed",
		RefStart:           0,
		OutputAsProportion: true,
		ResidualWindowSize: 1,
		Pseudocounts:       5,
		ExpectedFile:       "testdata/expected.myAln.highPseudo.bed",
		EnforceStrandMatch: false,
		GcContent:          0.5,
	},
	{MotifFile: "testdata/myMotifFile.txt",
		MotifType:          "Frequency",
		FastaFile:          "testdata/myAln.fa",
		ChromName:          "chr1",
		PropMatch:          0.1,
		OutFile:            "testdata/tmp.myAln.enforceStrand.bed",
		RefStart:           0,
		OutputAsProportion: true,
		ResidualWindowSize: 1,
		Pseudocounts:       0.1,
		ExpectedFile:       "testdata/expected.myAln.enforceStrand.bed",
		EnforceStrandMatch: true,
		GcContent:          0.5,
	},
	{MotifFile: "testdata/myMotifFile.txt",
		MotifType:          "Frequency",
		FastaFile:          "testdata/myAln.fa",
		ChromName:          "chr1",
		PropMatch:          0.8,
		OutFile:            "testdata/tmp.myAln.highPropMatch.bed",
		RefStart:           100,
		OutputAsProportion: true,
		ResidualWindowSize: 1,
		Pseudocounts:       0.1,
		ExpectedFile:       "testdata/expected.myAln.highPropMatch.bed",
		EnforceStrandMatch: false,
		GcContent:          0.5,
	},
	{MotifFile: "testdata/myMotifFile.txt",
		MotifType:          "Frequency",
		FastaFile:          "testdata/myAln.fa",
		ChromName:          "chr1",
		PropMatch:          0.1,
		OutFile:            "testdata/tmp.myAln.wideResidual.bed",
		RefStart:           100,
		OutputAsProportion: true,
		ResidualWindowSize: 50,
		Pseudocounts:       0.1,
		ExpectedFile:       "testdata/expected.myAln.wideResidual.bed",
		EnforceStrandMatch: false,
		GcContent:          0.5,
	},
}

func TestMatchComp(t *testing.T) {
	var err error
	var records []fasta.Fasta
	var s MatchCompSettings
	var epsilon float64 = 1e-6
	for _, v := range MatchCompTests {
		records = fasta.Read(v.FastaFile)
		fasta.AllToUpper(records)
		s = MatchCompSettings{
			MotifFile:          v.MotifFile,
			MotifType:          v.MotifType,
			Records:            records,
			PropMatch:          v.PropMatch,
			ChromName:          v.ChromName,
			OutFile:            v.OutFile,
			RefStart:           v.RefStart,
			Pseudocounts:       v.Pseudocounts,
			ResidualWindowSize: v.ResidualWindowSize,
			OutputAsProportion: v.OutputAsProportion,
			EnforceStrandMatch: v.EnforceStrandMatch,
			GcContent:          v.GcContent,
		}
		MatchComp(s)
		if !ApproxEquals(v.ExpectedFile, v.OutFile, epsilon) {
			t.Errorf("Error: Motif are not equal within a tolorance %v...", epsilon)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
