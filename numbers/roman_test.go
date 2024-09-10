package numbers

import (
	"testing"
)

func TestIntToRoman(t *testing.T) {
	testCases := []struct {
		input    int
		expected string
	}{
		{1, "I"},
		{3, "III"},
		{4, "IV"},
		{9, "IX"},
		{58, "LVIII"},
		{1994, "MCMXCIV"},
		{3994, "MMMCMXCIV"},
	}

	for _, test := range testCases {
		result := intToRoman(test.input)
		if result != test.expected {
			t.Errorf("Error: intToRoman(%d) = %s; expected = %s", test.input, result, test.expected)
		}
	}
}

func TestRomanToInt(t *testing.T) {
	testCases := []struct {
		input    string
		expected int
	}{
		{"I", 1},
		{"III", 3},
		{"IV", 4},
		{"IX", 9},
		{"LVIII", 58},
		{"MCMXCIII", 1993},
		{"MMMCMXCIV", 3994},
	}

	for _, test := range testCases {
		result := romanToInt(test.input)
		if result != test.expected {
			t.Errorf("Error: romanToInt(%s) = %d; expected = %d", test.input, result, test.expected)
		}
	}
}
