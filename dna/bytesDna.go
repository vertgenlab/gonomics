package dna

import (
	"github.com/vertgenlab/gonomics/common"
	"log"
	"strings"
)

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
	//VCF uses star to denote a deleted allele
	case '*':
		return Gap
	case '.':
		return Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", b)
		return N
	}
}

func ByteSliceToDnaBases(b []byte) []Base {
	var answer []Base = make([]Base, len(b))
	for i, byteValue := range b {
		answer[i] = ByteToBase(byteValue)
	}
	return answer
}

func ByteDnaBasesToString(dnaBases []Base) string {
	var str strings.Builder
	str.Grow(len(dnaBases))
	var err error
	for _, b := range dnaBases {
		_, err = str.WriteRune(BaseToRune(b))
		common.ExitIfError(err)
	}
	return str.String()
}
