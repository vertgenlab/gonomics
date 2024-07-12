package numbers

import (
	"math/rand"
	"testing"
	"time"
)

var seed *rand.Rand = rand.New(rand.NewSource(time.Now().UnixNano()))

func TestRandIntInRange(t *testing.T) {
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

func TestRandInt64InRange(t *testing.T) {
	x, y := int64(10), int64(20)
	for i := 0; i < 100; i++ {
		result := RandInt64InRange(x, y, seed)
		if result < x || result >= y {
			t.Errorf("RandInt64InRange(%d, %d) = %d; want a value in range [%d, %d)", x, y, result, x, y)
		}
	}
}

func TestRandFloat64InRange(t *testing.T) {
	x, y := float64(1.5), float64(3.5)
	for i := 0; i < 100; i++ {
		result := RandFloat64InRange(x, y, seed)
		if result < x || result >= y {
			t.Errorf("RandFloat64InRange(%f, %f) = %f; want a value in range [%f, %f)", x, y, result, x, y)
		}
	}
}
