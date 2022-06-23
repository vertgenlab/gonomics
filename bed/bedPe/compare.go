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
	if strings.Compare(a.Chrom1, b.Chrom1) != 0 {
		return false
	}
	if strings.Compare(a.Chrom2, b.Chrom2) != 0 {
		return false
	}
	if a.ChromStart1 != b.ChromStart1 {
		return false
	}
	if a.ChromStart2 != b.ChromStart2 {
		return false
	}
	if a.ChromEnd1 != b.ChromEnd1 {
		return false
	}
	if a.ChromEnd2 != b.ChromEnd2 {
		return false
	}
	return true
}
