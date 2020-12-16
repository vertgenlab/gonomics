package numbers

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"strings"
)

// EqualSliceInt is a function what will check two slices of int if the values are the same.
func EqualSliceInt(x []int, y []int) bool {
	if len(x) != len(y) {
		return false
	}
	for i := 0; i < len(x); i++ {
		if x[i] != y[i] {
			return false
		}
	}
	return true
}

// StringToInts will process a column of bytes and convert the slice into a slice of type int.
func StringToInts(column string) []int {
	work := strings.Split(column, ",")
	var answer []int = make([]int, len(work))
	for i := 0; i < len(work)-1; i++ {
		answer[i] = common.StringToInt(work[i])
	}
	return answer
}

// IntListToString will process a slice of type int as an input and return a each value separated by a comma as a string.
func IntListToString(nums []int) string {
	ans := strings.Builder{}
	ans.Grow(2 * len(nums))
	for i := 0; i < len(nums); i++ {
		ans.WriteString(fmt.Sprintf("%d", nums[i]))
		ans.WriteByte(',')
	}
	return ans.String()
}
