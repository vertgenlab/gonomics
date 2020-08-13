package numbers

import (
	"log"
	"math"
)

//these two values mark the overflow/underflow boundaries for math.Exp(x)
func LogCanConvert(x float64) bool {
	return x > 7.09782712893383973096e+02 || x < -7.45133219101941108420e+02
}

func AddLog(x float64, y float64) float64 {
	if math.IsInf(x, -1) {
		return y
	}
	if math.IsInf(y, -1) {
		return x
	}
	if x >= y {
		return x + math.Log(1+math.Exp(y-x))
	}
	return y + math.Log(1-math.Exp(y-x))
}

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
	return x + math.Log(1-math.Exp(y-x))
}

func MultiplyLog(x float64, y float64) float64 {
	if math.IsInf(x, -1) || math.IsInf(y, -1) {
		return math.Inf(-1)
	}
	return x + y
}

func DivideLog(x float64, y float64) float64 {
	if math.IsInf(x, -1) {
		return math.Inf(-1)
	}
	if math.IsInf(y, -1) {
		log.Fatalf("Divide by zero error in logSpace. x=%f. y=%f.", x, y)
	}
	return x - y
}
