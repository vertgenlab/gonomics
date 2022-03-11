package logspace

import (
	"log"
	"math"
)

// CanConvert returns true if the input logSpace number can be converted to a number in normal space without overflow or underflow.
func CanConvert(x float64) bool {
	return x < 709.4 && x > -745.1
}

// Average returns a log number that is the average of two input logSpace numbers
func Average(x float64, y float64) float64 {
	return Divide(Add(x, y), math.Log(2.0))
}

// Add returns the sum of two numbers in logSpace.
func Add(x float64, y float64) float64 {
	if math.IsInf(x, -1) {
		return y
	}
	if math.IsInf(y, -1) {
		return x
	}
	if x >= y {
		if CanConvert(y - x) {
			return x + math.Log1p(math.Exp(y-x))
		}
		return x
	}
	if CanConvert(y - x) {
		return y + math.Log1p(math.Exp(x-y))
	}
	return y
}

// Subtract returns the difference of two numbers in logSpace.
func Subtract(x float64, y float64) float64 {
	if x < y {
		log.Fatalf("Error: Taking the log of a negative number.")
		return 0
	}
	if x == y {
		return math.Inf(-1)
	}
	if math.IsInf(y, -1) {
		return x
	}
	if CanConvert(y - x) {
		return x + math.Log(1-math.Exp(y-x))
	}
	return x
}

// Multiply returns the product of two numbers in logSpace.
func Multiply(x float64, y float64) float64 {
	if math.IsInf(x, -1) || math.IsInf(y, -1) {
		return math.Inf(-1)
	}
	return x + y
}

// Divide returns the quotient of two numbers in logSpace.
func Divide(x float64, y float64) float64 {
	if math.IsInf(y, -1) {
		log.Fatalf("Divide by zero error in logSpace. x=%f. y=%f.", x, y)
	}
	if math.IsInf(x, -1) {
		return math.Inf(-1)
	}

	return x - y
}

// Pow returns log(exp(x)**y) where log is the natural logarithm.
// This function returns the log-space answer to x**y where x is already in log-space
// y is NOT in log-space
func Pow(x float64, y float64) float64 {
	// anything to the zero power is 1, so we return 0 == log(1).  0^0 is what this catches
	if y == 0.0 {
		return 0.0
	}
	return y * x
}
