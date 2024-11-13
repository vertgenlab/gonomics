package cigar

// equal asserts equality between two Cigar structs
func equal(a, b Cigar) bool {
	return a.Op == b.Op && a.RunLength == b.RunLength
}

// AllEqual performs a comparison between two slices of cigars
func AllEqual(a []Cigar, b []Cigar) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !equal(a[i], b[i]) {
			return false
		}
	}
	return true
}
