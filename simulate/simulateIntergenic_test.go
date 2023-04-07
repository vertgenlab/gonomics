package simulate

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"math/rand"
	"os"
	"testing"
)

var SimulateWithIndelsTests = []struct {
	FastaFile         string
	BranchLength      float64
	PropIndel         float64
	Lambda            float64
	GcContent         float64
	TransitionBias    float64
	VcfOutFile        string
	OutFastaFile      string
	ExpectedFastaFile string
	ExpectedVcfFile   string
	QName             string
}{
	{FastaFile: "testdata/rand.fa",
		BranchLength:      0.1,
		PropIndel:         0.2,
		Lambda:            1,
		GcContent:         0.42,
		TransitionBias:    1,
		VcfOutFile:        "testdata/tmp.vcf",
		OutFastaFile:      "testdata/tmp.rand.fa",
		ExpectedVcfFile:   "testdata/expected.rand.vcf",
		ExpectedFastaFile: "testdata/expected.rand.fa",
		QName:             "sim",
	},
	{FastaFile: "testdata/rand.fa",
		BranchLength:      0.1,
		PropIndel:         0.2,
		Lambda:            1,
		GcContent:         0.42,
		TransitionBias:    5,
		VcfOutFile:        "testdata/tmp.transition5.vcf",
		OutFastaFile:      "testdata/tmp.rand.transition5.fa",
		ExpectedVcfFile:   "testdata/expected.transition5.rand.vcf",
		ExpectedFastaFile: "testdata/expected.transition5.rand.fa",
		QName:             "sim",
	},
}

func TestSimulateWithIndels(t *testing.T) {
	rand.Seed(-1)
	var err error
	var records []fasta.Fasta
	for _, v := range SimulateWithIndelsTests {
		records = SimulateWithIndels(v.FastaFile, v.BranchLength, v.PropIndel, v.Lambda, v.GcContent, v.TransitionBias, v.VcfOutFile, v.QName)
		fasta.Write(v.OutFastaFile, records)
		if !fileio.AreEqual(v.OutFastaFile, v.ExpectedFastaFile) {
			t.Errorf("Error in SimulateWithIndels. Output fasta was not as expected.")
		} else {
			err = os.Remove(v.OutFastaFile)
			exception.PanicOnErr(err)
		}
		if v.VcfOutFile != "" {
			if !fileio.AreEqual(v.VcfOutFile, v.ExpectedVcfFile) {
				t.Errorf("Error in SimulateWithIndels. Output vcf was not as expected.")
			} else {
				err = os.Remove(v.VcfOutFile)
				exception.PanicOnErr(err)
			}
		}
	}
}
