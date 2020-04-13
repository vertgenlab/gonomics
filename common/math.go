package common

import (
	"log"
	"math/rand"
	"strconv"
)

func StringToInt(s string) int {
	n, err := strconv.Atoi(s)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a int\n", s)
	}
	return n
}

func StringToFloat64(s string) float64 {
	n, err := strconv.ParseFloat(s, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a float64\n", s)
	}
	return n
}

func StringToInt64(s string) int64 {
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to an int64\n", s)
	}
	return n
}

func StringToUint64(s string) uint64 {
	n, err := strconv.ParseUint(s, 10, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint64\n", s)
	}
	return n
}

func StringToUint32(s string) uint32 {
	n, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint32\n", s)
	}
	return uint32(n)
}

func Max(a int, b int) int {
	if a >= b {
		return a
	} else {
		return b
	}
}

func MaxInt32(a int32, b int32) int32 {
	if a >= b {
		return a
	} else {
		return b
	}
}

func MaxUint32(a uint32, b uint32) uint32 {
	if a >= b {
		return a
	} else {
		return b
	}
}

func MaxInt64(a int64, b int64) int64 {
	if a >= b {
		return a
	} else {
		return b
	}
}

func MaxFloat64(a float64, b float64) float64 {
	if a >= b {
		return a
	} else {
		return b
	}
}

func Min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}

func MinInt32(a int32, b int32) int32 {
	if a <= b {
		return a
	} else {
		return b
	}
}

func MinUint32(a uint32, b uint32) uint32 {
	if a <= b {
		return a
	} else {
		return b
	}
}

func MinInt64(a int64, b int64) int64 {
	if a <= b {
		return a
	} else {
		return b
	}
}

func MinFloat64(a float64, b float64) float64 {
	if a <= b {
		return a
	} else {
		return b
	}
}

func TripleMax(a int, b int, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= c {
		return b
	} else {
		return c
	}
}

func TripleMin(a int, b int, c int) int {
	if a <= b && a <= c {
		return a
	} else if b <= c {
		return b
	} else {
		return c
	}
}

func RandIntInRange(x int, y int) int {
	return int(rand.Float64()*float64(y-x)) + x
}

func RandInt64InRange(x int64, y int64) int64 {
	return int64(rand.Float64()*float64(y-x)) + x
}
