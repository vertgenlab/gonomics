package dna

import (
	"bytes"
	"fmt"
	"log"
	"unicode/utf8"
)

func runeToBase(r rune) (Base, error) {
	switch r {
	case 'A':
		return A, nil
	case 'C':
		return C, nil
	case 'G':
		return G, nil
	case 'T':
		return T, nil
	case 'N':
		return N, nil
	case 'a':
		return a, nil
	case 'c':
		return c, nil
	case 'g':
		return g, nil
	case 't':
		return t, nil
	case 'n':
		return n, nil
	case '-':
		return Gap, nil
	default:
		return N, fmt.Errorf("dna: unexpected character in dna %r", r)
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
	default:
		log.Fatal(fmt.Errorf("dna: unexpected value in dna Base when converting to rune %u", base))
		return 'N'
	}
}

func StringToBases(s string) ([]Base, error) {
	answer := make([]Base, utf8.RuneCountInString(s))
	var err error

	for index, runeValue := range s {
		answer[index], err = runeToBase(runeValue)
		if err != nil {
			return nil, err
		}
	}
	return answer, nil
}

func BasesToString(bases []Base) string {
	var buffer bytes.Buffer

	for _, base := range bases {
		buffer.WriteRune(BaseToRune(base))
	}
	return buffer.String()
}
