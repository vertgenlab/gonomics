package pFasta

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"testing"
)

var WriteTests = []struct {
	OutFile      string
	Records      []pFasta
	ExpectedFile string
}{
	{OutFile: "testdata/out.test.pFa",
		Records: []pFasta{
			pFasta{Name: "chr1",
				Seq: []pDna.Float32Base{
					pDna.Float32Base{
						A: 0.3,
						C: 0.5,
						G: 0.2,
						T: 0,
					},
					pDna.Float32Base{
						A: 0,
						C: 0.99,
						G: 0.01,
						T: 0,
					},
				},
			},
			pFasta{Name: "chr2",
				Seq: []pDna.Float32Base{
					pDna.Float32Base{ //gap base
						A: 0,
						C: 0,
						G: 0,
						T: 0,
					},
					pDna.Float32Base{
						A: 0,
						C: 0,
						G: 0.2,
						T: 0.8,
					},
					pDna.Float32Base{
						A: 0.25,
						C: 0.25,
						G: 0.25,
						T: 0.25,
					},
				},
			},
		},
		ExpectedFile: "testdata/expected.test.pFa",
	},
}

func TestWriteAndRead(t *testing.T) {
	var records []pFasta
	for _, v := range WriteTests {
		Write(v.OutFile, v.Records)
		records = Read(v.OutFile)
		fmt.Println(records)
	}
}
