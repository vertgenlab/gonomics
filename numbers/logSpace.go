package numbers

import (
	"log"
	"math"
	//DEBUG: "fmt"
)

//LogCanConvert returns true if the input logSpace number can be converted to a number in normal space without overflow or underflow.
func LogCanConvert(x float64) bool {
	return x < 709.4 && x > -745.1
}

//AverageLog returns the average of a list of logSpace numbers
func AverageLog(x []float64) float64 {
	var n int = len(x)
	var sum float64 = x[0]
	var logN float64 = math.Log(float64(n))
	for i := 1; i < len(x); i++ {
		sum = AddLog(sum, x[i])
	}
	return DivideLog(sum, logN)
}

//MidpointLog returns a log number that is equidistant from two input logSpace numbers
func MidpointLog(x float64, y float64) float64 {
	return DivideLog(AddLog(x, y), math.Log(2.0))
}

//AddLog returns the sum of two numbers in logSpace.
func AddLog(x float64, y float64) float64 {
	if math.IsInf(x, -1) {
		return y
	}
	if math.IsInf(y, -1) {
		return x
	}
	if x >= y {
		if LogCanConvert(y - x) {
			return x + math.Log1p(math.Exp(y-x))
		}
		return x
	}
	if LogCanConvert(y - x) {
		return y + math.Log1p(math.Exp(x-y))
	}
	return y
}

//SubtractLog returns the difference of two numbers in logSpace.
func SubtractLog(x float64, y float64) float64 {
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
	if LogCanConvert(y - x) {
		return x + math.Log(1-math.Exp(y-x))
	}
	return x
}

//MultiplyLog returns the product of two numbers in logSpace.
func MultiplyLog(x float64, y float64) float64 {
	if math.IsInf(x, -1) || math.IsInf(y, -1) {
		return math.Inf(-1)
	}
	return x + y
}

//DivideLog returns the quotient of two numbers in logSpace.
func DivideLog(x float64, y float64) float64 {
	if math.IsInf(y, -1) {
		log.Fatalf("Divide by zero error in logSpace. x=%f. y=%f.", x, y)
	}
	if math.IsInf(x, -1) {
		return math.Inf(-1)
	}

	return x - y
}

//LogPow returns log(x**y) where log is the natural logarithm. Safe for large numbers. Support for positive real numbers with integer exponents.
func LogPow(x float64, y float64) float64 {
	if x < 0 {
		log.Fatalf("LowPowInt does not handle negative x values. x=%e\n", x)
	}
	// anything to the zero power is 1, so we return 0 == log(1).  0^0 is what this catches
	if y == 0.0 {
		return 0.0
	}
	return y * math.Log(x)
}

// PowLog returns log(exp(x)**y) where log is the natural logarithm.
// This back the log-space answer to x**y where x is already in log-space
// In other words, this function returns the log-space answer to x**y where x is already in log-space
func PowLog(x float64, y float64) float64 {
	// anything to the zero power is 1, so we return 0 == log(1).  0^0 is what this catches
	if y == 0.0 {
		return 0.0
	}
	return y * x
}
