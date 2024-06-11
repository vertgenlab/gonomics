package numbers

import (
	"math/rand"
)

// RandIntInRange returns a pseudorandom value of type int between x and y.
// Output includes x, but not y.
func RandIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandIntInRangeSrc(x int, y int, src rand.Source) int {
	return int(src.Int63())%(y-x+1) + x
}

// RandInt64InRange returns a pseudorandom value of type int64 between x and y.
// Output includes x, but not y.
func RandInt64InRange(x int64, y int64) int64 {
	return int64(rand.Float64()*float64(y-x)) + x
}

// RandFloat64InRange returns a pseudorandom value of type int between x and y.
func RandFloat64InRange(x float64, y float64) float64 {
	return rand.Float64()*(y-x) + x
}
