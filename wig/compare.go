package wig

import (
	"strings"
)

//isEqual returns true if two Wig data structures contain the exact same data and returns false otherwise.
func isEqual(alpha *Wig, beta *Wig) bool {
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
			if alpha.Values[i].Position != beta.Values[i].Position {
				return false
			}
			if alpha.Values[i].Value != beta.Values[i].Value {
				return false
			}
		}
	}
	return true
}

//AllEqual returns true if two input arrays of wig data structures contain all the same data, false otherwise.
func AllEqual(alpha []*Wig, beta []*Wig) bool {
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
