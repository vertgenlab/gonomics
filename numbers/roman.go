package numbers

import (
	"strings"
)

// intToRoman converts an integer to its Roman numeral representation returning an integer.
func intToRoman(numeral int) string {
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

// romanToInt converts a Roman numeral string to its corresponding integer value.
func romanToInt(roman string) int {
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
