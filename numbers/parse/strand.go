package parse

import (
	"log"
)

// StringToStrand returns true if s is '+', false if '-', and fatal error otherwise.
func StringToStrand(s string) bool {
	switch s {
	case "+":
		return true
	case "-":
		return false
	default:
		log.Fatalf("Error: expecting %s to be a strand that is either '+' or '-'\n", s)
		return false
	}
}

// StrandToRune returns '+' if true, and otherwise '-'.
func StrandToRune(strand bool) rune {
	if strand {
		return '+'
	} else {
		return '-'
	}
}
