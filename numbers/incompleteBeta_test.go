package numbers

import (
	"testing"
)

var incompleteBetaTests = []struct {
	a      float64
	b      float64
	x      float64
	answer float64
}{
	{0.1, 0.1, 0.5, 0.5},
	{0.5, 0.5, 0.5, 0.5},
	{0.1, 0.1, 0.2, 0.4397092},
	{10, 10, 0.1, 3.929882e-06},
	{10, 10, 0.3, 0.03255336},
	{10, 10, 0.5, 0.50000000},
	{10, 10, 0.7, 0.9674466},
	{10, 10, 1.0, 1.0000000},
	{15, 10, 0.5, 0.1537281},
	{15, 10, 0.6, 0.4890802},
	{10, 15, 0.5, 0.8462719},
	{10, 15, 0.6, 0.9783419},
	{20, 20, 0.4, 0.1020586},
	{40, 40, 0.4, 0.03581031},
	{40, 40, 0.7, 0.9999017},
}

func TestIncompleteBeta(t *testing.T) {
	for _, test := range incompleteBetaTests {
		calculated := incompleteBetaHelper(test.a, test.b, test.x)
		if calculated > test.answer*1.01 || calculated < test.answer*0.99 {
			t.Errorf("For incompleteBetaHelper with a: %f, b: %f, and x: %f we expected %e, but got %e.", test.a, test.b, test.x, test.answer, calculated)
		}
	}
}
