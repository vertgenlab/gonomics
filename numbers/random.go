package numbers

import (
	"math/rand"
)

func RandIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandInt64InRange(x int64, y int64) int64 {
	return int64(rand.Float64()*float64(y-x)) + x
}

func RandFloat64InRange(x float64, y float64) float64 {
	return rand.Float64()*(y-x) + x
}
