package numbers

import (
	"math"
)

// roundSigFigs rounds a float64 to a specified number of significant figures and returns as a float32
func RoundSigFigs(num float64, sigFigs int) float32 {
	if (num == 0) {
		return float32(0)
	}
	if (num < 0) {num = -num}
	answer := math.Ceil(math.Log10(num))
	power := sigFigs - int(answer)
	magnitude := math.Pow10(power)
	shift := math.Round(num*magnitude)
	res := float32(shift/magnitude)
	return res
}