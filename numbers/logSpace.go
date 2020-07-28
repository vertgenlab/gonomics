package numbers

import (
	"math"
	"log"
)

//these two values mark the overflow/underflow boundaries for math.Exp(x)
func LogCanConvert(x float64) bool {
	return x > 7.09782712893383973096e+02 || x < -7.45133219101941108420e+02
}

func AddLog(x float64, y float64) float64 {
	//if x == -infinity {
	//	return y
	//}
	//if y == -infinity {
	//	return x
	//}
	if x >= y {
		return x + math.Log(1 + math.Exp(y-x))
	}
	return y + math.Log(1 - math.Exp(y-x))
}

func SubtractLog(x float64, y float64) float64 {
	if x < y {
		log.Fatalf("Error: Taking the log of a negative number.")
		return 0
	} 
	/*if x == y{
		return -Inf
	} 
	if y == -Inf {
		return x
	}*/
	return x + math.Log(1 - math.Exp(y-x))
}

func MultiplyLog(x float64, y float64) float64 {
	/*if x == -Inf || y == -Inf {
		return -Inf
	}*/
	return x + y
}