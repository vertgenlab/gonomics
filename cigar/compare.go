package cigar

// isEqual checks if two Cigar structs are equal.
func isEqual(a Cigar, b Cigar) bool {
	return (a.RunLength == b.RunLength && a.Op == b.Op)
}

// isByteCigarEqual checks if two ByteCigar structs are equal.
func isByteCigarEqual(a ByteCigar, b ByteCigar) bool {
	return (a.RunLen == b.RunLen && a.Op == b.Op)
}

// EqualByteCigar checks if two slices of ByteCigar are equal.
func EqualByteCigar(a []ByteCigar, b []ByteCigar) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !isByteCigarEqual(a[i], b[i]) {
			return false
		}
	}
	return true
}
