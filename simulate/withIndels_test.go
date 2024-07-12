package simulate

import (
	"math/rand"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var randFaInput string = "testdata/rand.fa.gz"

var WithIndelsTests = []struct {
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
	{FastaFile: randFaInput,
		BranchLength:      0.1,
		PropIndel:         0.2,
		Lambda:            1,
		GcContent:         0.42,
		TransitionBias:    1,
		VcfOutFile:        "testdata/tmp.vcf.gz",
		ExpectedVcfFile:   "testdata/expected.rand.vcf.gz",
		ExpectedFastaFile: "testdata/expected.rand.fa.gz",
		QName:             "sim",
	},
	{FastaFile: randFaInput,
		BranchLength:      0.1,
		PropIndel:         0.2,
		Lambda:            1,
		GcContent:         0.42,
		TransitionBias:    5,
		VcfOutFile:        "testdata/tmp.transition5.vcf.gz",
		ExpectedVcfFile:   "testdata/expected.transition5.rand.vcf.gz",
		ExpectedFastaFile: "testdata/expected.transition5.rand.fa.gz",
		QName:             "sim",
	},
}

func TestWithIndels(t *testing.T) {
	seed := rand.New(rand.NewSource(-1))
	var records []fasta.Fasta
	for _, v := range WithIndelsTests {
		records = WithIndels(v.FastaFile, v.BranchLength, v.PropIndel, v.Lambda, v.GcContent, v.TransitionBias, v.VcfOutFile, v.QName, seed)
		if !fasta.AllAreEqual(records, fasta.Read(v.ExpectedFastaFile)) {
			fasta.Write(v.ExpectedFastaFile, records)
			t.Errorf("Error in SimulateWithIndels. Output fasta was not as expected.")
		}

		if v.VcfOutFile != "" {
			if !fileio.AreEqual(v.VcfOutFile, v.ExpectedVcfFile) {
				exception.PanicOnErr(os.Rename(v.VcfOutFile, v.ExpectedVcfFile))
				t.Errorf("Error in SimulateWithIndels. Output vcf was not as expected.")
			} else {
				fileio.EasyRemove(v.VcfOutFile)
			}
		}
	}
}
