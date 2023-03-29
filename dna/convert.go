package dna

import (
	"errors"
	"fmt"
	"log"
	"strings"
	"unicode"

	"github.com/vertgenlab/gonomics/exception"
)

// RuneToBase converts a rune into a dna.Base if it matches one of the acceptable DNA characters.
// Note: '*', used by VCF to denote deleted alleles becomes Nil
func RuneToBase(r rune) (Base, error) {
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
		return LowerA, nil
	case 'c':
		return LowerC, nil
	case 'g':
		return LowerG, nil
	case 't':
		return LowerT, nil
	case 'n':
		return LowerN, nil
	case '-':
		return Gap, nil
	case '*':
		return Nil, nil
	case '.':
		return Dot, nil
	default:
		return N, errors.New(fmt.Sprintf("ERROR: '%s' is an invalid base. This error is caused by invalid characters in a DNA sequence. "+
			"Note that gonomics does not support extended IUPAC nucleotide codes (only AaCcGgTtNn-). If you input a fasta file, "+
			"check to make sure the sequences only use valid characters. Invalid characters can be masked with Ns using the faFormat command in gonomics.", string(r)))
	}
}

// ByteToBase converts a byte into a dna.Base if it matches one of the acceptable DNA characters.
// Notes: It will also mask the lower case values and return dna.Base as uppercase bases.
// Note: '*', used by VCF to denote deleted alleles, becomes a Gap in DNA.
func ByteToBase(b byte) (Base, error) {
	switch b {
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
		return LowerA, nil
	case 'c':
		return LowerC, nil
	case 'g':
		return LowerG, nil
	case 't':
		return LowerT, nil
	case 'n':
		return LowerN, nil
	case '-':
		return Gap, nil
	case '*':
		return Nil, nil
	case '.':
		return Dot, nil
	default:
		return N, errors.New(fmt.Sprintf("ERROR: '%s' is an invalid base. This error is caused by invalid characters in a DNA sequence. "+
			"Note that gonomics does not support extended IUPAC nucleotide codes (only AaCcGgTtNn-). If you input a fasta file, "+
			"check to make sure the sequences only use valid characters. Invalid characters can be masked with Ns using the faFormat command in gonomics.", string(b)))
	}
}

// BaseToRune converts a dna.Base type into a rune.
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
		log.Panicf("Error: unexpected value in dna Base when converting to rune\n")
		return 'N'
	}
}

// Extract returns a subsequence of an input slice of DNA bases from an input start and end point.
func Extract(rec []Base, start int, end int) []Base {
	// Use the 3-arg slicing for safety: any append will force a copy instead of overwriting rec
	return rec[start:end:end]
}

// BaseToString converts a DNA base to a string by casting a BaseToRune result to a string.
func BaseToString(b Base) string {
	return string(baseToByteArray[b])
}

// StringToBase parses a string into a single DNA base
func StringToBase(s string) Base {
	if len(s) != 1 {
		log.Panic("String with a length other than 1 can not be turned into a dna.Base\n")
	}
	var b Base
	var err error
	b, err = ByteToBase(s[0])
	exception.PanicOnErr(err)
	return b
}

// StringToBases parses a string into a slice of DNA bases
func StringToBases(s string) []Base {
	answer := make([]Base, len(s))
	var b Base
	var err error
	for index := range s {
		b, err = ByteToBase(s[index])
		exception.PanicOnErr(err)
		answer[index] = b
	}
	return answer
}

// StringToBasesForced parses a string into a slice of DNA bases and N-masks any invalid characters.
func StringToBasesForced(s string) []Base {
	answer := make([]Base, len(s))
	var b Base
	var err error
	for index := range s {
		b, err = ByteToBase(s[index])
		if err != nil {
			if unicode.IsUpper(rune(s[index])) {
				b = N
			} else {
				b = LowerN
			}
		}
		answer[index] = b
	}
	return answer
}

// baseToByte is an efficient lookup for the rune corresponding to a given dna.Base.
// intended to remain as a private array to help the BasesToString function.
// panics if value input is not a valid Base.
// quicker than BaseToByte by ~5x
var baseToByteArray = []byte{'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', '-', '.', '*'}

// BasesToString converts a slice of DNA bases into a string. Useful for writing to files.
func BasesToString(bases []Base) string {
	var buffer strings.Builder
	buffer.Grow(len(bases))
	var err error
	for i := range bases {
		err = buffer.WriteByte(baseToByteArray[bases[i]])
		if err != nil {
			log.Panicf("problem writing byte")
		}
	}
	return buffer.String()
}

// ByteSliceToDnaBases will convert a slice of bytes into a slice of Bases.
func ByteSliceToDnaBases(b []byte) []Base {
	var answer []Base = make([]Base, len(b))
	var d Base
	var err error
	for i := range b {
		d, err = ByteToBase(b[i])
		exception.PanicOnErr(err)
		answer[i] = d
	}
	return answer
}
