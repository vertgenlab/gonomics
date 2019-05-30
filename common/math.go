package common

import (
	"fmt"
	"strconv"
)

func StringToInt64(s string) int64 {
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		ExitIfError(fmt.Errorf("Error: trouble converting %s to an int64\n", s))
	}
	return n
}

func Max(a int, b int) int {
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
