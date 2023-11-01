package bedpe

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
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

// AnnotateFeetDist calculates the distance between two bedpe feet and puts the value in the annotation field
func AnnotateFeetDist(b []BedPe) {
	var dist int
	for i := range b {
		b[i].A.FieldsInitialized = 11
		dist = numbers.Max(b[i].A.ChromStart, b[i].B.ChromStart) - numbers.Min(b[i].A.ChromStart, b[i].B.ChromStart)
		b[i].A.Annotation = append(b[i].A.Annotation, fileio.IntToString(dist))
	}
}
