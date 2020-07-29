package dnaThreeBit

import (
	"bytes"
	"log"
	"unicode/utf8"
)

func RuneToThreeBitBase(r rune) ThreeBitBase {
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
		return A
	case 'c':
		return C
	case 'g':
		return G
	case 't':
		return T
	case 'n':
		return N
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", r)
		return N
	}
}

func ThreeBitBaseToRune(base Base) rune {
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
	default:
		log.Fatalf("Error: unexpected value in dna Base when converting to rune\n")
		return 'N'
	}
}

func ThreeBitBaseToString(b ThreeBitBase) string {
	return string(ThreeBitBaseToRune(b))
}

func StringToBases(s string) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index, runeValue := range s {
		answer[index] = RuneToBase(runeValue)
	}
	return answer
}

func ToString(bases []Base) string {
	var buffer bytes.Buffer

	for i:=0; i<fragment.Len; i++ {
		buffer.WriteRune(dna.BaseToRune(GetBase(fragment, i)))
	}
	return buffer.String()
}

func ToDnaBases(fragment *ThreeBit) []dna.Base {
        answer := make([]dna.Base, fragment.Len)
        for i:=0; i<fragment.Len; i++ {
                answer[i] = GetBase(fragment, i)
        }
        return answer
}

