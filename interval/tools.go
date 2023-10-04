package interval

import (
	"fmt"
	"github.com/vertgenlab/gonomics/numbers"
)

// CoordsToString takes in a interval and returns a string in the format of chr:start-end
func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}

func IntervalSize(i Interval) int {
	return i.GetChromEnd() - i.GetChromStart()
}

func OverlapSize(a, b Interval) int {
	if a.GetChrom() != b.GetChrom() {
		return 0
	}
	end := numbers.Min(a.GetChromEnd(), b.GetChromEnd())
	start := numbers.Max(a.GetChromStart(), b.GetChromStart())
	// overlapSize := end - start
	if overlapSize <= 0 {
		return 0
	} else {
		return overlapSize
	}
}
