package bedpe

import (
	"strings"
)

func AllAreEqual(a []BedPe, b []BedPe) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

// Equal returns true if two input BedPe entries have the same coordinates for both regions. False otherwise.
func Equal(a BedPe, b BedPe) bool {
	if strings.Compare(a.ChromA, b.ChromB) != 0 {
		return false
	}
	if strings.Compare(a.ChromB, b.ChromB) != 0 {
		return false
	}
	if a.ChromStartA != b.ChromStartA {
		return false
	}
	if a.ChromStartB != b.ChromStartB {
		return false
	}
	if a.ChromEndA != b.ChromEndA {
		return false
	}
	if a.ChromEndB != b.ChromEndB {
		return false
	}
	return true
}
