package interval

import (
	"fmt"
)

// CoordsToString takes in a interval and returns a string in the format of chr:start-end.
func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}
