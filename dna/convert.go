package dna

import (
	"strings"
)

/*
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
		return LowerA
	case 'c':
		return LowerC
	case 'g':
		return LowerG
	case 't':
		return LowerT
	case 'n':
		return LowerN
	case '-':
		return Gap
	//VCF uses star to denote a deleted allele
	case '*':
		return Gap
	case '.':
		return Dot
	default:
		log.Fatalf("Error: unexpected character in dna %v\n", r)
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
		log.Fatalf("Error: unexpected character in dna %v\n", b)
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
	case LowerA:
		return 'a'
	case LowerC:
		return 'c'
	case LowerG:
		return 'g'
	case LowerT:
		return 't'
	case LowerN:
		return 'n'
	case Gap:
		return '-'
	case Dot:
		return '.'
	default:
		log.Fatalf("Error: unexpected value in dna Base when converting to rune: %s\n", string(base))
		return 'N'
	}
}

//Extract returns a subsequence of an input slice of DNA bases from an input start and end point.
func Extract(rec []Base, start int64, end int64) []Base {
	return rec[start:end]
}
*/

//BaseToString converts a DNA base to a string by casting a BaseToRune result to a string.
func BaseToString(b Base) string {
	return string(b)
}

//StringToBases parses a string into a slice of DNA bases
func StringToBases(s string) []Base {
	return []Base(s)
}

//BasesToString converts a slice of DNA bases into a string. Useful for writing to files.
func BasesToString(bases []Base) string {
	var s strings.Builder
	s.Grow(len(bases))
	for i := range bases {
		s.WriteByte(byte(bases[i]))
	}
	return s.String()
}

// ByteSliceToDnaBases will convert a slice of bytes into a slice of Bases with no lowercase bases.
func ByteSliceToDnaBases(b []byte) []Base {
	answer := make([]Base, len(b))
	for i := range b {
		answer[i] = Base(b[i])
	}
	return answer
}
