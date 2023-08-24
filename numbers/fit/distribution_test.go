package fit

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"testing"
)

var ZeroInflatedNegativeBinomialTests = []struct {
	R             float64
	P             float64
	Z             float64
	NumDataPoints int
	LearningRate  float64
	Delta         float64
	Epsilon       float64
}{
	{
		R:             2,
		P:             0.8,
		Z:             0,
		NumDataPoints: 10_000,
		LearningRate:  1e-6,
		Delta:         1e-6,
		Epsilon:       1e-6,
	},
}

func TestZeroInflatedNegativeBinomial(t *testing.T) {
	var currR, currP, currZ float64
	for _, v := range ZeroInflatedNegativeBinomialTests {
		data := simulateZeroInflatedNegativeBinomial(v.R, v.P, v.Z, v.NumDataPoints)
		fmt.Println(data)
		currR, currP, currZ = ZeroInflatedNegativeBinomial(data, v.LearningRate, v.Delta, v.Epsilon)
		if currR != v.R {
			t.Errorf("Error: currR: %v. expectedR: %v.\n", currR, v.R)
		}
		if currP != v.P {
			t.Errorf("Error: currP: %v. ExpectedP: %v.\n", currP, v.P)
		}
		if currZ != v.Z {
			t.Errorf("Error: currZ: %v. ExpectedZ: %v.\n", currZ, v.Z)
		}
	}
}

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
	var lines []string
	var currR float64
	var currP float64
	var couldNotEvaluate bool
	for _, v := range NegativeBinomialTests {
		variates = make([]float64, 0)
		lines = fileio.Read(v.InFile)
		for i := range lines {
			variates = append(variates, parse.StringToFloat64(lines[i]))
		}
		currR, currP, couldNotEvaluate = NegativeBinomial(variates)
		if couldNotEvaluate {
			log.Fatalf("Error: Could not evaluate distribution for given values.\n")
		}
		if (currR-v.ExpectedR)/v.ExpectedR > 0.05 {
			t.Errorf("Error: currR: %v not as expected: %v\n", currR, v.ExpectedR)
		}
		if (currP-v.ExpectedP)/v.ExpectedP > 0.05 {
			t.Errorf("Error: currP: %v not as expected: %v\n", currP, v.ExpectedP)
		}
	}
}

var PoissonTests = []struct {
	InFile         string
	ExpectedLambda float64
}{
	{InFile: "testdata/poissonVariates.lambda.1.txt.gz",
		ExpectedLambda: 1,
	},
	{InFile: "testdata/poissonVariates.lambda.2.txt.gz",
		ExpectedLambda: 2,
	},
	{InFile: "testdata/poissonVariates.lambda.3.txt.gz",
		ExpectedLambda: 3,
	},
	{InFile: "testdata/poissonVariates.lambda.4.txt.gz",
		ExpectedLambda: 4,
	},
	{InFile: "testdata/poissonVariates.lambda.5.txt.gz",
		ExpectedLambda: 5,
	},
}

func TestPoisson(t *testing.T) {
	var variates []float64
	var lines []string
	var currLambda float64
	for _, v := range PoissonTests {
		variates = make([]float64, 0)
		lines = fileio.Read(v.InFile)
		for i := range lines {
			variates = append(variates, parse.StringToFloat64(lines[i]))
		}
		currLambda = Poisson(variates)
		if (currLambda-v.ExpectedLambda)/v.ExpectedLambda > 0.05 {
			t.Errorf("Error: currLambda: %v is not as expected: %v\n", currLambda, v.ExpectedLambda)
		}
	}
}
