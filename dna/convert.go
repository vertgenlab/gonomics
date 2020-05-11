package dna

import (
	"bytes"
	"log"
	"unicode/utf8"
)

func RuneToBase(r rune) Base {
	switch r {
	case 'A':
		return A
	case 'C':
		return C
	case 'G':
		return G
	case 'T':
		return T
	case 'N':
		return N
	case 'a':
		return a
	case 'c':
		return c
	case 'g':
		return g
	case 't':
		return t
	case 'n':
		return n
	case '-':
		return Gap
	//VCF uses star to denote a deleted allele
	case '*':
		return Gap
	case '.':
		return Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", r)
		return N
	}
}

func BaseToRune(base Base) rune {
	switch base {
	case A:
		return 'A'
	case C:
		return 'C'
	case G:
		return 'G'
	case T:
		return 'T'
	case N:
		return 'N'
	case a:
		return 'a'
	case c:
		return 'c'
	case g:
		return 'g'
	case t:
		return 't'
	case n:
		return 'n'
	case Gap:
		return '-'
	case Dot:
		return '.'
	default:
		log.Fatalf("Error: unexpected value in dna Base when converting to rune\n")
		return 'N'
	}
}

func Extract(rec []Base, start int64, end int64) []Base {
	return rec[start:end]
}

func BaseToString(b Base) string {
	return string(BaseToRune(b))
}

func StringToBases(s string) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index, runeValue := range s {
		answer[index] = RuneToBase(runeValue)
	}
	return answer
}

func BasesToString(bases []Base) string {
	var buffer bytes.Buffer

	for _, base := range bases {
		buffer.WriteRune(BaseToRune(base))
	}
	return buffer.String()
}
