package genePred

import (
	"strings"
)

func AllAreEqual(a []GenePred, b []GenePred) (bool, []int) {
	var fields []int
	var equal bool

	if len(a) != len(b) {
		return false, fields
	}
	for i := 0; i < len(a); i++ {
		equal, fields = Equal(a[i], b[i])
		return equal, fields
	}
	return equal, fields
}

func Equal(a GenePred, b GenePred) (bool, []int) {
	var isEqual = true
	var fields = make([]int, 9)

	if strings.Compare(a.Id, b.Id) != 0 {
		fields[0] = 1
		isEqual = false
	}
	if strings.Compare(a.Chrom, b.Chrom) != 0 {
		fields[1] = 1
		isEqual = false
	}
	if a.Strand != b.Strand {
		fields[2] = 1
		isEqual = false
	}
	if a.TxStart != b.TxStart {
		fields[3] = 1
		isEqual = false
	}
	if a.TxEnd != b.TxEnd {
		fields[4] = 1
		isEqual = false
	}
	if a.CdsStart != b.CdsStart {
		fields[5] = 1
		isEqual = false
	}
	if a.CdsEnd != b.CdsEnd {
		fields[6] = 1
		isEqual = false
	}
	//exon ends must have the same number of values as exon starts
	for i := 0; i < len(a.ExonStarts); i++ {
		if a.ExonStarts[i] != b.ExonStarts[i] {
			fields[8] = 1
			isEqual = false
		}
		if b.ExonEnds[i] != b.ExonEnds[i] {
			fields[9] = 1
			isEqual = false
		}
	}
	return isEqual, fields
}
