package pFasta

import (
	"math/rand"
	"testing"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fileio"
)

func randSeq(length int) PFasta {
	var one, two, three, four, total float32
	var answer PFasta = PFasta{}
	answer.Name = "randSeq"
	answer.Seq = make([]pDna.Float32Base, length)
	for i := 0; i < length; i++ {
		one = rand.Float32()
		two = rand.Float32()
		three = rand.Float32()
		four = rand.Float32()
		total = one + two + three + four
		answer.Seq[i].A = one / total
		answer.Seq[i].C = two / total
		answer.Seq[i].G = three / total
		answer.Seq[i].T = four / total
	}
	return answer
}

var WriteTests = []struct {
	OutFile   string
	Records   []PFasta
	Precision float32
}{
	{OutFile: "testdata/out.test.pFa",
		Records: []PFasta{
			{Name: "chr1",
				Seq: []pDna.Float32Base{
					{
						A: 0.3,
						C: 0.5,
						G: 0.2,
						T: 0,
					},
					{
						A: 0,
						C: 0.99,
						G: 0.01,
						T: 0,
					},
				},
			},
			{Name: "chr2",
				Seq: []pDna.Float32Base{
					{ //gap base
						A: 0,
						C: 0,
						G: 0,
						T: 0,
					},
					{
						A: 0,
						C: 0,
						G: 0.2,
						T: 0.8,
					},
					{
						A: 0.25,
						C: 0.25,
						G: 0.25,
						T: 0.25,
					},
				},
			},
		},
		Precision: 1e-3, // float16 for writing is quite inaccurate, this fails at 1e-4
	},
	{
		OutFile:   "testdata/out.test.pFa",
		Records:   []PFasta{randSeq(100000)},
		Precision: 1e-2,
	},
}

func TestWriteAndRead(t *testing.T) {
	var records []PFasta
	for _, v := range WriteTests {
		Write(v.OutFile, v.Records)
		records = Read(v.OutFile)
		if !AllAreEqual(records, v.Records, v.Precision) {
			t.Errorf("Error: in pFasta. Write and read test was not as expected.\n")
		} else {
			fileio.MustRemove(v.OutFile)
		}
	}
}
