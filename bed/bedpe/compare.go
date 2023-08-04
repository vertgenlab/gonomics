package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"sort"
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
	if !bed.Equal(a.A, b.A) {
		return false
	}
	if !bed.Equal(a.B, b.B) {
		return false
	}
	return true
}

// SortByCoord will return a bedpe sorted by the A half
func SortByCoord(bedpeFile []BedPe) {
	sort.Slice(bedpeFile, func(i, j int) bool { return Compare(bedpeFile[i], bedpeFile[j]) == -1 })
}

// Compare returns zero for equal BedPes and otherwise returns the ordering of the two BedPe entries. Used for SortByCoord.
func Compare(one BedPe, two BedPe) int {
	chromComp := strings.Compare(one.A.Chrom, two.A.Chrom)
	if chromComp != 0 {
		return chromComp
	}
	if one.A.ChromStart < two.A.ChromStart {
		return -1
	}
	if one.A.ChromStart > two.A.ChromStart {
		return 1
	}
	if one.A.ChromEnd < two.A.ChromEnd {
		return -1
	}
	if one.A.ChromEnd > two.A.ChromEnd {
		return 1
	}
	return 0
}
