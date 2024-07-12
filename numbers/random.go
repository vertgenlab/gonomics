package numbers

import (
	"math/rand"
)

// RandIntInRange returns a pseudorandom value of type int between x and y.
// Output includes x, but not y.
func RandIntInRange(x int, y int, seed *rand.Rand) int {
	return seed.Intn(y-x) + x
}

// RandInt64InRange returns a pseudorandom value of type int64 between x and y.
// Output includes x, but not y.
func RandInt64InRange(x int64, y int64, seed *rand.Rand) int64 {
	return int64(seed.Float64()*float64(y-x)) + x
}

// RandFloat64InRange returns a pseudorandom value of type int between x and y.
func RandFloat64InRange(x float64, y float64, seed *rand.Rand) float64 {
	return seed.Float64()*(y-x) + x
}
