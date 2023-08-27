package interval

import (
	"fmt"
)

// CoordsToString takes in a interval and returns a string in the format of chr:start-end
func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}

func AreEqual(a Interval, b Interval) bool {
	var isEqual bool = false

	if a.GetChrom() == b.GetChrom() {
		if a.GetChromStart() == b.GetChromStart() {
			if a.GetChromEnd() == b.GetChromEnd() {
				isEqual = true
			}
		}
	}

	return isEqual
}
