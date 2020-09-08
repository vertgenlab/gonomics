package dna

import (
	"github.com/vertgenlab/gonomics/common"
	"log"
	"strings"
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

// ByteToBaseToBase converts a byte into a dna.Base if it matches one of the acceptable DNA characters.
// Notes: It will also mask the lower case values and return dna.Base as uppercase bases.
// Note: '*', used by VCF to denote deleted alleles, becomes a Gap in DNA.
func ByteToBase(b byte) Base {
	switch b {
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
		return A
	case 'c':
		return C
	case 'g':
		return G
	case 't':
		return T
	case 'n':
		return N
	case '-':
		return Gap
	case '*':
		return Gap
	case '.':
		return Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", b)
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
	var buffer strings.Builder
	buffer.Grow(len(bases))
	var err error
	for _, b := range bases {
		_, err = buffer.WriteRune(BaseToRune(b))
		common.ExitIfError(err)
	}
	return buffer.String()
}

// ByteSliceToDnaBases will convert a slice of bytes into a slice of Bases with no lowercase bases.
func ByteSliceToDnaBases(b []byte) []Base {
	var answer []Base = make([]Base, len(b))
	for i, byteValue := range b {
		answer[i] = ByteToBase(byteValue)
	}
	return answer
}
