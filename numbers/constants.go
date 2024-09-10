package numbers

import (
	"math/bits"
)

const MaxInt = 1<<(bits.UintSize-1) - 1
const MinInt = -MaxInt - 1
const MaxUint = 1<<bits.UintSize - 1

var numerals = []int{1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1}
var romans = []string{"M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I"}

var symbols = map[byte]int{
	'I': 1,
	'V': 5,
	'X': 10,
	'L': 50,
	'C': 100,
	'D': 500,
	'M': 1000,
}
