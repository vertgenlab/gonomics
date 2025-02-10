package parse

import (
	"log"
	"strconv"
)

// StringToBool parses a string into a bool and will exit on error.
func StringToBool(s string) bool {
	b, err := strconv.ParseBool(s)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a bool", s)
	}
	return b
}

// StringToInt parses a string into an int and will exit on error.
func StringToInt(s string) int {
	// because strings for ints can be written in scientific notation,
	// instead of converting from string directly to int (cannot handle scientific notation),
	// convert string to float64 (nFloatIntermediate), then to int
	nFloatIntermediate, err := strconv.ParseFloat(s, 64)
	n := int(nFloatIntermediate)
	// to check that the string (e.g. "3" or "3.5") was indeed an int (e.g. 3) and not a non-int number (e.g. 3.5),
	// after the string (e.g. "3.5") to float64 (e.g. nFloatIntermediate == 3.5) to int (e.g. 3) conversion,
	// convert int to float64 (e.g. nFloatIntermediate2 == 3) again
	// if nFloatIntermediate2 (3) != nFloatIntermediate (3.5), then the string was a non-int number, error
	nFloatIntermediate2 := float64(n)
	if nFloatIntermediate2 != nFloatIntermediate {
		log.Panicf("Error: StringToInt was used, but the String \"%s\" was not Int", s)
	}
	//n, err := strconv.Atoi(s) // previous code for converting string directly to int
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a int", s)
	}
	return n
}

// StringToFloat32 parses a string into a float32 and will exit on error.
func StringToFloat32(s string) float32 {
	n, err := strconv.ParseFloat(s, 32)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a float32\n", s)
	}
	return float32(n)
}

// StringToFloat64 parses a string into a float64 and will exit on error.
func StringToFloat64(s string) float64 {
	n, err := strconv.ParseFloat(s, 64)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a float64\n", s)
	}
	return n
}

// StringToInt8 parses a string into an int8 and will exit on error.
func StringToInt8(s string) int8 {
	n, err := strconv.ParseInt(s, 10, 8)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a int8\n", s)
	}
	return int8(n)
}

// StringToInt16 parses a string into an int16 and will exit on error.
func StringToInt16(s string) int16 {
	n, err := strconv.ParseInt(s, 10, 16)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a int16\n", s)
	}
	return int16(n)
}

// StringToInt32 parses a string into an int32 and will exit on error.
func StringToInt32(s string) int32 {
	n, err := strconv.ParseInt(s, 10, 32)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a int32\n", s)
	}
	return int32(n)
}

// StringToInt64 parses a string into an int64 and will exit on error.
func StringToInt64(s string) int64 {
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a int64\n", s)
	}
	return n
}

// StringToUint64 parses a string into a uint64 and will exit on error.
func StringToUint64(s string) uint64 {
	n, err := strconv.ParseUint(s, 10, 64)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a uint64\n", s)
	}
	return n
}

// StringToUint32 parses a string into a uint32 and will exit on error.
func StringToUint32(s string) uint32 {
	n, err := strconv.ParseUint(s, 10, 32)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a uint32\n", s)
	}
	return uint32(n)
}

// StringToUint16 parses a string into a uint16 and will exit on error.
func StringToUint16(s string) uint16 {
	n, err := strconv.ParseUint(s, 10, 16)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a uint16\n", s)
	}
	return uint16(n)
}

// StringToUint8 parses a string into a uint8 and will exit on error.
func StringToUint8(s string) uint8 {
	n, err := strconv.ParseUint(s, 10, 8)
	if err != nil {
		log.Panicf("Error: trouble converting \"%s\" to a uint8\n", s)
	}
	return uint8(n)
}
