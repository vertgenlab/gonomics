package bedpe

import (
	"github.com/vertgenlab/gonomics/numbers"
	"testing"
)

func TestFindStats(t *testing.T) {
	var x float64
	norm := []float64{3.524983097071441e-06, 4.8341057860438124e-05, 6.9173998139866e-05}
	math := Read("testdata/mathTest.bedpe")
	SortByCoord(math)
	a, mu, s := FindStats(math)
	for i := range a {
		x = numbers.NormalDist(a[i], mu[i], s[i])
		if x != norm[i] {
			t.Errorf("Normal Distribution calculation does not match expected value")
		}

	}

}
