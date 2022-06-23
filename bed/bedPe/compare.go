package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
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
	if !bed.Equal(a.A, b.A) {
		return false
	}
	if !bed.Equal(a.B, b.B) {
		return false
	}
	return true
}
