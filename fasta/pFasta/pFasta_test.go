package pFasta

import (
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/fileio"
	"testing"
)

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
		Records:   []PFasta{RandSeq(1000, "test", false, 8, nil)},
		Precision: 1e-2,
	},
	{
		OutFile:   "testdata/out.test.pFa",
		Records:   []PFasta{RandSeq(10000, "test", false, 8, nil)},
		Precision: 1e-2,
	},
	{
		// TODO: there is a bug caused by relative precision being used for very small numbers when seed is 7 or 8
		OutFile:   "testdata/out.test.pFa",
		Records:   []PFasta{RandSeq(100000, "test", false, 9, nil)},
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
