package wig

import (
	"sort"
	"strings"
)

//isEqual returns true if two Wig data structures contain the exact same data and returns false otherwise.
func isEqual(alpha Wig, beta Wig) bool {
	if strings.Compare(alpha.StepType, beta.StepType) != 0 {
		return false
	}
	if strings.Compare(alpha.Chrom, beta.Chrom) != 0 {
		return false
	}
	if alpha.Start != beta.Start {
		return false
	}
	if alpha.Step != beta.Step {
		return false
	}
	if len(alpha.Values) == len(beta.Values) {
		for i := 0; i < len(alpha.Values); i++ {
			if alpha.Values[i] != beta.Values[i] {
				return false
			}
		}
	}
	return true
}

//AllEqual returns true if two input arrays of wig data structures contain all the same data, false otherwise.
func AllEqual(alpha []Wig, beta []Wig) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			return false
		}
	}
	return true
}

//Compare returns zero for equal wigs and otherwise returns the ordering of the two wig entries. Used for SortByCoord.
func Compare(alpha Wig, beta Wig) int {
	chromComp := strings.Compare(alpha.Chrom, beta.Chrom)
	if chromComp != 0 {
		return chromComp
	}
	if alpha.Start < beta.Start {
		return -1
	}
	if alpha.Start > beta.Start {
		return 1
	}
	if len(alpha.Values) < len(beta.Values) {
		return -1
	}
	if len(alpha.Values) > len(beta.Values) {
		return 1
	}
	return 0
}

//SortByCoord sorts in place a slice of Wig structs by their genomic position.
func SortByCoord(w []Wig) {
	sort.Slice(w, func(i, j int) bool { return Compare(w[i], w[j]) == -1 })
}
