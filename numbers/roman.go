package numbers

import (
	"strings"
)

// IntToRoman converts an integer to its Roman numeral representation returning an integer.
func IntToRoman(numeral int) string {
	var answer strings.Builder
	for numeral > 0 {
		for i := range numerals {
			if numeral >= numerals[i] {
				answer.WriteString(romans[i])
				numeral -= numerals[i]
				break
			}
		}
	}
	return answer.String()
}

// RomanToInt converts a Roman numeral string to its corresponding integer value.
func RomanToInt(roman string) int {
	var answer, prev, curr int
	for i := len(roman) - 1; i >= 0; i-- {
		curr = symbols[roman[i]]
		if curr < prev {
			answer -= curr
		} else {
			answer += curr
		}
		prev = curr
	}
	return answer
}
