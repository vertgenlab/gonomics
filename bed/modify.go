package bed

import "log"

//Trim shortens bed entries on the left and right side by an input-specified number of bases. These values must not exceed the length of the bed entry and must be non-negative.
func Trim(b []Bed, trimLeft int, trimRight int) {
	for i := range b {
		if trimLeft < 0 || trimRight < 0 {
			log.Fatalf("Error in bed/Trim. Must trim bed values by a value greater or equal to zero.")
		}
		b[i].ChromStart = b[i].ChromStart + trimLeft
		b[i].ChromEnd = b[i].ChromEnd - trimRight
		if b[i].ChromStart >= b[i].ChromEnd {
			log.Fatalf("Error in Trim. Attempted to remove too much from bed entry. Please select a lower trim value or exclude the bed entry as position %v\t%v.\n", b[i].Chrom, b[i].ChromStart)
		}
	}
}
