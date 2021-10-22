package bed

import "log"

//Trim
func Trim(b []Bed, trimLeft int, trimRight int) {
	for i := range b {
		b[i].ChromStart = b[i].ChromStart + trimLeft
		b[i].ChromEnd = b[i].ChromEnd - trimRight
		if b[i].ChromStart >= b[i].ChromEnd {
			log.Fatalf("Error in Trim. Attempted to remove too much from bed entry. Please select a lower trim value or exclude the bed entry as position %v\t%v.\n", b[i].Chrom, b[i].ChromStart)
		}
	}
}
