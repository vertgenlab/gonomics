package common

import (
	"log"
)

// true if '+', false if '-', and fatal error otherwise
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

func StrandToRune(strand bool) rune {
	if strand {
		return '+'
	} else {
		return '-'
	}
}
