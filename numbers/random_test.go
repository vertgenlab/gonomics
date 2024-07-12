package numbers

import (
	"math/rand"
	"testing"
	"time"
)

var seed *rand.Rand = rand.New(rand.NewSource(time.Now().UnixNano()))

func TestCompareGlobalLocalSource(t *testing.T) {
	// Initialize the global random source
	rand.New(rand.NewSource(time.Now().UnixNano()))
	global := rand.Intn(1000)

	// Initialize the local random source
	local := RandIntInRange(0, 1000, rand.New(rand.NewSource(time.Now().UnixNano())))

	// Ensure the global and local random sources produce different outputs
	if global == local {
		t.Errorf("Error: Global and local random sources produced the same value: %v", global)
	}
}

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
				t.Errorf("Error: Value %d is outside the expected range [%d, %d]", v, interval.x, interval.y)
			}
		}
	}
}

func TestRandInt64InRange(t *testing.T) {
	x, y := int64(10), int64(20)
	for i := 0; i < 10; i++ {
		result := RandInt64InRange(x, y, seed)
		if result < x || result >= y {
			t.Errorf("Error: RandInt64InRange(%d, %d) = %d; want a value in range [%d, %d)", x, y, result, x, y)
		}
	}
}

func TestRandFloat64InRange(t *testing.T) {
	x, y := float64(1.5), float64(3.5)
	for i := 0; i < 10; i++ {
		result := RandFloat64InRange(x, y, seed)
		if result < x || result >= y {
			t.Errorf("Error: RandFloat64InRange(%f, %f) = %f; want a value in range [%f, %f)", x, y, result, x, y)
		}
	}
}
