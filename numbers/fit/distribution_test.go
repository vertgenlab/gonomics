package fit

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"testing"
)

var NegativeBinomialTests = []struct {
	InFile    string
	ExpectedR float64
	ExpectedP float64
}{
	{InFile: "testdata/negativeBinomialVariates.r.1.p.0.05.txt.gz",
		ExpectedR: 1,
		ExpectedP: 0.05,
	},
	{InFile: "testdata/negativeBinomialVariates.r.1.p.0.1.txt.gz",
		ExpectedR: 1,
		ExpectedP: 0.10,
	},
	{InFile: "testdata/negativeBinomialVariates.r.1.p.0.15.txt.gz",
		ExpectedR: 1,
		ExpectedP: 0.15,
	},
	{InFile: "testdata/negativeBinomialVariates.r.1.p.0.2.txt.gz",
		ExpectedR: 1,
		ExpectedP: 0.2,
	},
	{InFile: "testdata/negativeBinomialVariates.r.1.p.0.25.txt.gz",
		ExpectedR: 1,
		ExpectedP: 0.25,
	},
	{InFile: "testdata/negativeBinomialVariates.r.2.p.0.05.txt.gz",
		ExpectedR: 2,
		ExpectedP: 0.05,
	},
	{InFile: "testdata/negativeBinomialVariates.r.2.p.0.1.txt.gz",
		ExpectedR: 2,
		ExpectedP: 0.10,
	},
	{InFile: "testdata/negativeBinomialVariates.r.2.p.0.15.txt.gz",
		ExpectedR: 2,
		ExpectedP: 0.15,
	},
	{InFile: "testdata/negativeBinomialVariates.r.2.p.0.2.txt.gz",
		ExpectedR: 2,
		ExpectedP: 0.2,
	},
	{InFile: "testdata/negativeBinomialVariates.r.2.p.0.25.txt.gz",
		ExpectedR: 2,
		ExpectedP: 0.25,
	},
	{InFile: "testdata/negativeBinomialVariates.r.4.p.0.05.txt.gz",
		ExpectedR: 4,
		ExpectedP: 0.05,
	},
	{InFile: "testdata/negativeBinomialVariates.r.4.p.0.1.txt.gz",
		ExpectedR: 4,
		ExpectedP: 0.10,
	},
	{InFile: "testdata/negativeBinomialVariates.r.4.p.0.15.txt.gz",
		ExpectedR: 4,
		ExpectedP: 0.15,
	},
	{InFile: "testdata/negativeBinomialVariates.r.4.p.0.2.txt.gz",
		ExpectedR: 4,
		ExpectedP: 0.2,
	},
	{InFile: "testdata/negativeBinomialVariates.r.4.p.0.25.txt.gz",
		ExpectedR: 4,
		ExpectedP: 0.25,
	},
}

func TestNegativeBinomial(t *testing.T) {
	var variates []float64
	var currR float64
	var currP float64
	for _, v := range NegativeBinomialTests {
		variates = make([]float64, 0)
		lines := fileio.Read(v.InFile)
		for i := range lines {
			variates = append(variates, parse.StringToFloat64(lines[i]))
			currR, currP = NegativeBinomial(variates)
		}
		if (currR-v.ExpectedR)/v.ExpectedR > 0.05 {
			t.Errorf("Error: currR: %v not as expected: %v\n", currR, v.ExpectedR)
		}
		if (currP-v.ExpectedP)/v.ExpectedP > 0.05 {
			t.Errorf("Error: currP: %v not as expected: %v\n", currP, v.ExpectedP)
		}
	}
}
