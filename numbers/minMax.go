package numbers

//Max returns the maximum of two int inputs of type int.
func Max(a int, b int) int {
	if a >= b {
		return a
	} else {
		return b
	}
}

//MaxInt32 returns the maximum of two inputs of type int32.
func MaxInt32(a int32, b int32) int32 {
	if a >= b {
		return a
	} else {
		return b
	}
}

//MaxUint32 returns the maximum of two inputs of type uint32.
func MaxUint32(a uint32, b uint32) uint32 {
	if a >= b {
		return a
	} else {
		return b
	}
}

//MaxInt64 returns the maximum of two inputs of type int64.
func MaxInt64(a int64, b int64) int64 {
	if a >= b {
		return a
	} else {
		return b
	}
}

//MaxFloat64 returns the maximum of two inputs of type float64.
func MaxFloat64(a float64, b float64) float64 {
	if a >= b {
		return a
	} else {
		return b
	}
}

//Min returns the minimum of tw inputs of type int.
func Min(a int, b int) int {
	if a <= b {
		return a
	} else {
		return b
	}
}

//MinInt32 returns the minimum of two inputs of type int32.
func MinInt32(a int32, b int32) int32 {
	if a <= b {
		return a
	} else {
		return b
	}
}

//MinUint32 returns the minimum of two inputs of type uint32.
func MinUint32(a uint32, b uint32) uint32 {
	if a <= b {
		return a
	} else {
		return b
	}
}

//MinInt64 returns the minimum of two inputs of type int64.
func MinInt64(a int64, b int64) int64 {
	if a <= b {
		return a
	} else {
		return b
	}
}

//MinFloat64 returns the minimum of two inputs of type float64.
func MinFloat64(a float64, b float64) float64 {
	if a <= b {
		return a
	} else {
		return b
	}
}

//TripleMax returns the maximum of three inputs of type int.
func TripleMax(a int, b int, c int) int {
	if a >= b && a >= c {
		return a
	} else if b >= c {
		return b
	} else {
		return c
	}
}

//TripeMin returns the minimum of three inputs of type int.
func TripleMin(a int, b int, c int) int {
	if a <= b && a <= c {
		return a
	} else if b <= c {
		return b
	} else {
		return c
	}
}

//MaxIntSlice returns the maximum of a slice of type int.
func MaxIntSlice(a []int) int {
	var answer int = a[0]
	for i := 1; i < len(a); i++ {
		if a[i] > answer {
			answer = a[i]
		}
	}
	return answer
}

//MinIntSlice returns the minimum of a slice of type int.
func MinIntSlice(a []int) int {
	var answer int = a[0]
	for i := 1; i < len(a); i++ {
		if a[i] < answer {
			answer = a[i]
		}
	}
	return answer
}
