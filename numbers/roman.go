package numbers

import "strings"

// intToRoman converts an integer to its Roman numeral representation.
func intToRoman(num int) string {
	var builder strings.Builder
	for num > 0 {
		for i := range numerals {
			if num >= numerals[i] {
				builder.WriteString(romans[i])
				num -= numerals[i]
				break
			}
		}
	}
	return builder.String()
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
