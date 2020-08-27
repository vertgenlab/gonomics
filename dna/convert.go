package dna

import (
	"bytes"
	"log"
	"unicode/utf8"
)

//RuneToBase converts a rune into a dna.Base if it matches one of the acceptable DNA characters.
//Note: '*', used by VCF to denote deleted alleles, becomes a Gap in DNA.
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

//BaseToRune converts a dna.Base type into a rune.
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

//Extract returns a subsequence of an input slice of DNA bases from an input start and end point.
func Extract(rec []Base, start int64, end int64) []Base {
	return rec[start:end]
}

//BaseToString converts a DNA base to a string by casting a BaseToRune result to a string.
func BaseToString(b Base) string {
	return string(BaseToRune(b))
}

//StringToBases parses a string into a slice of DNA bases
func StringToBases(s string) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index, runeValue := range s {
		answer[index] = RuneToBase(runeValue)
	}
	return answer
}

//BasesToString converts a slice of DNA bases into a string. Useful for writing to files.
func BasesToString(bases []Base) string {
	var buffer bytes.Buffer

	for _, base := range bases {
		buffer.WriteRune(BaseToRune(base))
	}
	return buffer.String()
}
