package simulate

import (
	"fmt"
	"math/rand"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

var GCcontent = 0.42

var RandGeneTests = []struct {
	name   string
	length int
	GC     float64
}{
	{"testingRandGene", 102, GCcontent},
}

func TestRandGene(t *testing.T) {
	rand.Seed(1)
	for _, test := range RandGeneTests {
		a := RandGene(test.name, test.length, test.GC)
		if len(a[0].Seq) != test.length {
			t.Errorf("expected RandGene to give %v, gave %v", test.length, len(a[0].Seq))
		}
		//fmt.Print(dna.BasesToString(a[0].Seq), "\n")
		//fasta.Write("outputRandGene.fasta", a)
	}
}

var MutateSeqTests = []struct {
	sequence     string
	branchLength float64
	gp           string
}{
	{"testdata/longDebug.fasta", 0.5, "testdata/longDebug.gp"}, //branch length of 1 gives higher chance of returning a new base so you can see a difference even with a short sequence
}

func TestMutateGene(t *testing.T) {
	rand.Seed(1)
	for _, test := range MutateSeqTests {
		seq := fasta.Read(test.sequence)
		bases := seq[0].Seq
		a := MutateGene(bases, test.branchLength, test.gp, true)
		if len(bases) != len(a) {
			t.Errorf("Expected same length sequences. Original: %v \n Ending: %v", len(bases), len(a))
		}
	}
}

var IndelLengthTests = []struct {
	Lambda       float64
	OutFile      string
	ExpectedFile string
}{
	{1, "testdata/test.IndelLength.Lambda1.txt", "testdata/expected.IndelLength.Lambda1.txt"},
	{0.5, "testdata/test.IndelLength.LambdaPoint5.txt", "testdata/expected.IndelLength.LambdaPoint5.txt"},
	{3, "testdata/test.IndelLength.Lambda3.txt", "testdata/expected.IndelLength.Lambda3.txt"},
}

func TestIndelLength(t *testing.T) {
	var variateCount = 10000
	var err error
	rand.Seed(23)
	var lengths = make([]string, variateCount)
	var i int

	for _, v := range IndelLengthTests {
		i = 0
		for i < variateCount {
			lengths[i] = fmt.Sprintf("%v", indelLength(v.Lambda))
			i++
		}
		fileio.Write(v.OutFile, lengths)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Stability error in IndelLengths.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
