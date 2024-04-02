package wig

import (
	"log"
	//DEBUG: "fmt".
	"math/rand"
	"sort"
	"strings"

	"github.com/vertgenlab/gonomics/numbers"
)

// isEqual returns true if two Wig data structures contain the exact same data values and returns false otherwise.
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

// AllEqual returns true if two slices of Wig data structures contain all the same data vales in the same order, false otherwise.
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

// compare returns zero for equal wigs and otherwise returns the ordering of the two wig entries.
// Order is based on start position and then by length of values. Used for SortByCoord.
func compare(alpha Wig, beta Wig) int {
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

// SortByCoord sorts in place a slice of Wig structs by their genomic position with secondary sorting by the number of values.
func SortByCoord(w []Wig) {
	sort.Slice(w, func(i, j int) bool { return compare(w[i], w[j]) == -1 })
}

// Pearson calculates the Pearson Correlation Coefficient between two input wig slices. Data values equal to the 'missing' value are ignored.
// samplingFrequency is a number between 0 and 1. This represents the proportion of bases that are considered for evaluating the PCC.
func Pearson(alpha map[string]Wig, beta map[string]Wig, missing float64, samplingFrequency float64) float64 {
	var a = make([]float64, 0)
	var b = make([]float64, 0)
	var currPos int
	var currKey string
	var foundInMap bool
	var r float64
	if samplingFrequency < 0 || samplingFrequency > 1 {
		log.Fatalf("Error in wig.Pearson. samplingFrequency must be a value between 0 and 1.")
	}
	for currKey = range alpha {
		if _, foundInMap = beta[currKey]; !foundInMap {
			log.Fatalf("Error: in wig Pearson, chromin the first wig: %s, not found in the second wig.\n", currKey)
		}
		if len(alpha[currKey].Values) != len(beta[currKey].Values) {
			log.Fatalf("Error in wig.Pearson. Entries with the same chr name have a different number of values. Len(a): %v. Len(b): %v.", len(alpha[currKey].Values), len(beta[currKey].Values))
		}
		for currPos = range alpha[currKey].Values {
			if alpha[currKey].Values[currPos] != missing && beta[currKey].Values[currPos] != missing { //if neither wig has the missing value at this position.
				r = rand.Float64()
				if r < samplingFrequency {
					a = append(a, alpha[currKey].Values[currPos])
					b = append(b, beta[currKey].Values[currPos])
				}
			}
		}
	}
	return numbers.Pearson(a, b)
}
