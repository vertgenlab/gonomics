package popgen

import (
	"github.com/vertgenlab/gonomics/exception"
	"testing"
)

var MaximumLikelihoodTests = []struct {
	VcfFile                 string
	Left                    float64
	Right                   float64
	Error                   float64
	UnPolarized             bool
	DivergenceAscertainment bool
	D                       int
	IntegralError           float64
	Verbose                 int
	ExpectedValue           float64
}{
	{VcfFile: "testdata/simulated.alpha4.N100.S100.seed19.vcf",
		Left:                    -10,
		Right:                   10,
		Error:                   1e-5,
		UnPolarized:             true,
		DivergenceAscertainment: false,
		IntegralError:           1e-5,
		Verbose:                 0,
		ExpectedValue:           3.893977236832579,
	},
}

func TestMaximumLikelihood(t *testing.T) {
	var testValue float64
	var data *Afs
	var err error
	var s MleSettings
	for _, v := range MaximumLikelihoodTests {
		s = MleSettings{
			Left:                    v.Left,
			Right:                   v.Right,
			Error:                   v.Error,
			UnPolarized:             v.UnPolarized,
			DivergenceAscertainment: v.DivergenceAscertainment,
			D:                       v.D,
			IntegralError:           v.IntegralError,
			Verbose:                 v.Verbose,
		}
		data, err = VcfToAfs(v.VcfFile, v.UnPolarized, v.DivergenceAscertainment, false)
		exception.PanicOnErr(err)
		testValue = SelectionMaximumLikelihoodEstimate(*data, s)
		if testValue != v.ExpectedValue {
			t.Errorf("Error in MaximumLikelihood. Value %v is not as expected.", testValue)
		}
	}
}
