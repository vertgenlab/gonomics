package numbers

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
