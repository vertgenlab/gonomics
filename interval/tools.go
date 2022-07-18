package interval

import (
	"fmt"
)

func CoordsToString(i Interval) string {
	return fmt.Sprintf("%s:%v-%v", i.GetChrom(), i.GetChromStart(), i.GetChromEnd())
}
