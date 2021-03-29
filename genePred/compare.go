package genePred

import (
	"strings"
)

// AllAreEqual examines if two genePred files are the same entry by entry.
func AllAreEqual(a []GenePred, b []GenePred) (bool, []int) {
	var fields []int
	var equal bool

	if len(a) != len(b) {
		return false, fields
	}
	for i := 0; i < len(a); i++ {
		equal, fields = Equal(a[i], b[i])
		if !equal {
			return equal, fields
		}
	}
	return equal, fields
}

// Equal checks if two genePred entires are identical. This will return a bool of whether the whole entries are identical, and a slice of ints that refer to which fields did not match.
//Ex: if one or more exonStart positions do not match fields will contain a 9.
func Equal(a GenePred, b GenePred) (bool, []int) {
	var isEqual = true
	var fields []int

	if strings.Compare(a.Id, b.Id) != 0 {
		fields = append(fields, 1)
		isEqual = false
	}
	if strings.Compare(a.Chrom, b.Chrom) != 0 {
		fields = append(fields, 2)
		isEqual = false
	}
	if a.Strand != b.Strand {
		fields = append(fields, 3)
		isEqual = false
	}
	if a.TxStart != b.TxStart {
		fields = append(fields, 4)
		isEqual = false
	}
	if a.TxEnd != b.TxEnd {
		fields = append(fields, 5)
		isEqual = false
	}
	if a.CdsStart != b.CdsStart {
		fields = append(fields, 6)
		isEqual = false
	}
	if a.CdsEnd != b.CdsEnd {
		fields = append(fields, 7)
		isEqual = false
	}
	//exon ends must have the same number of values as exon starts
	for i := 0; i < len(a.ExonStarts); i++ {
		if a.ExonStarts[i] != b.ExonStarts[i] {
			fields = append(fields, 9)
			isEqual = false
		}
		if a.ExonEnds[i] != b.ExonEnds[i] {
			fields = append(fields, 10)
			isEqual = false
		}
	}
	return isEqual, fields
}
