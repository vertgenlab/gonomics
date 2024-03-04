package numbers

import (
	"testing"
)

var roundToSigFigsTests = []struct {
	Num		float64
	SigFigs	int
	Expected float32	
}{
	{Num: 0.2837562,
	SigFigs: 4,
	Expected: 0.2838,
		},
	{Num: 0.273942,
		SigFigs: 7,
		Expected: 0.273942,
		},
	{Num: 1.2837523,
		SigFigs: 3,
		Expected: 1.28,
		},	
}

func TestRoundToSigFigs(t *testing.T) {
	for _, testCase := range roundToSigFigsTests {
		res := RoundSigFigs(testCase.Num, testCase.SigFigs)
		if res != testCase.Expected {
			t.Errorf("Error: in in pFasta. roundToSigFigs test not as expected.")
		}
	}
}
