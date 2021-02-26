package numbers

import (
	"math/bits"
)

const MaxInt = 1<<(bits.UintSize-1) - 1
const MinInt = -MaxInt - 1
const MaxUint = 1<<bits.UintSize - 1
