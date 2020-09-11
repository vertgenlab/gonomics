package common

import (
	"log"
	"math/rand"
	"strconv"
)

// StringToInt parses a string into an int and will exit on error
func StringToInt(s string) int {
	n, err := strconv.Atoi(s)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a int\n", s)
	}
	return n
}

// StringToFloat64 parses a string into a float64 and will exit on error
func StringToFloat64(s string) float64 {
	n, err := strconv.ParseFloat(s, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a float64\n", s)
	}
	return n
}

// StringToInt8 parses a string into an int8 and will exit on error
func StringToInt8(s string) int8 {
	n, err := strconv.ParseInt(s, 10, 8)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a int8\n", s)
	}
	return int8(n)
}

// StringToInt16 parses a string into an int16 and will exit on error
func StringToInt16(s string) int16 {
	n, err := strconv.ParseUint(s, 10, 16)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a int16\n", s)
	}
	return int16(n)
}

// StringToInt32 parses a string into an int32 and will exit on error
func StringToInt32(s string) int32 {
	n, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a int32\n", s)
	}
	return int32(n)
}

// StringToInt64 parses a string into an int64 and will exit on error
func StringToInt64(s string) int64 {
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to an int64\n", s)
	}
	return n
}

// StringToUint64 parses a string into a uint64 and will exit on error
func StringToUint64(s string) uint64 {
	n, err := strconv.ParseUint(s, 10, 64)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint64\n", s)
	}
	return n
}

// StringToUint32 parses a string into a uint32 and will exit on error
func StringToUint32(s string) uint32 {
	n, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint32\n", s)
	}
	return uint32(n)
}

// StringToUint16 parses a string into a uint16 and will exit on error
func StringToUint16(s string) uint16 {
	n, err := strconv.ParseUint(s, 10, 16)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint16\n", s)
	}
	return uint16(n)
}

// StringToUint8 parses a string into a uint8 and will exit on error
func StringToUint8(s string) uint8 {
	n, err := strconv.ParseUint(s, 10, 8)
	if err != nil {
		log.Fatalf("Error: trouble converting %s to a uint8\n", s)
	}
	return uint8(n)
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
