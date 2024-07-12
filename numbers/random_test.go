package numbers

import (
	"math/rand"
	"testing"
)

func TestRandIntInRange(t *testing.T) {
	seed := rand.New(rand.NewSource(0))
	n := 20
	values := make([]int, n)
	tests := []struct {
		x, y int
	}{
		{0, 1},
		{0, 10},
		{-5, 5},
		{100, 200},
	}
	for _, interval := range tests {
		for i := 0; i < n; i++ {
			values[i] = RandIntInRange(interval.x, interval.y, seed)
		}
		for _, v := range values {
			if v < interval.x || v > interval.y {
				t.Errorf("Value %d is outside the expected range [%d, %d]", v, interval.x, interval.y)
			}
		}
	}
}
