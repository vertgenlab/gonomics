package numbers

import (
	"math/rand"
	"testing"
	"time"
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

func TestRandInt64InRange(t *testing.T) {
	seed := rand.New(rand.NewSource(time.Now().UnixNano()))
	n := 100
	values := make([]int64, n)
	tests := []struct {
		x, y int64
	}{
		{0, 1},
		{0, 99},
		{-5, 5},
		{-1, 1},
	}
	for _, interval := range tests {
		for i := 0; i < n; i++ {
			values[i] = RandInt64InRange(interval.x, interval.y, seed)
		}
		for _, v := range values {
			if v < interval.x || v > interval.y {
				t.Errorf("Value %d is outside the expected range [%d, %d]", v, interval.x, interval.y)
			}
		}
	}
}

func TestRandFloat64InRange(t *testing.T) {
	seed := rand.New(rand.NewSource(0))
	var samples int = 30
	values := make([]float64, samples)
	tests := []struct {
		name string
		x, y float64
	}{
		{
			x: 10.0,
			y: 20.0,
		},
		{
			x: -5.0,
			y: 5.0,
		},
		{
			x: 3.14159,
			y: 203.14159,
		},
	}
	for _, interval := range tests {
		for i := range values {
			values[i] = RandFloat64InRange(interval.x, interval.y, seed)
		}
		// Check if all values are within the range
		for _, v := range values {
			if v <= interval.x || v >= interval.y { // Note: '>' changed to '>=' for float comparison
				t.Errorf("Value %.4f is outside the expected range [%.4f, %.4f]", v, interval.x, interval.y)
			}
		}
	}
}
