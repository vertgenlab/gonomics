package dnaThreeBit

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
)

// RuneToThreeBitBase returns a single bases in ThreeBitBase format that corresponds to the given rune
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

// ThreeBitBaseToRune returns a rune that corresponds to the single base in ThreeBitBase format
func ThreeBitBaseToRune(base ThreeBitBase) rune {
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

// ThreeBitBaseToString returns a string that corresponds to the single base give as a ThreeBitBase
func ThreeBitBaseToString(b ThreeBitBase) string {
	return string(ThreeBitBaseToRune(b))
}

// FromString creates a new ThreeBit from a string of DNA characters {A,C,G,T,N}
func FromString(s string) *ThreeBit {
	answer := &ThreeBit{Seq: []uint64{}, Len: 0}
	for _, runeValue := range s {
		answer = Append(answer, RuneToThreeBitBase(runeValue))
	}
	return answer
}

// ToString returns a string representation of the ThreeBit passed in
func ToString(fragment *ThreeBit) string {
	var buffer strings.Builder
	buffer.Grow(fragment.Len)
	for i := 0; i < fragment.Len; i++ {
		buffer.WriteRune(dna.BaseToRune(GetBase(fragment, i)))
	}
	return buffer.String()
}

// RangeToDnaBases returns a slice of dna.Base that represents the bases from
// start to end (left-closed, right-open) of fragment
func RangeToDnaBases(fragment *ThreeBit, start int, end int) []dna.Base {
	if end > fragment.Len || start >= end {
		log.Fatalf("Error: unable to extract bases from %d to %d from a sequence of length %d\n", start, end, fragment.Len)
	}
	answer := make([]dna.Base, 0, end-start)
	for i := start; i < end; i++ {
		answer = append(answer, GetBase(fragment, i))
	}
	return answer
}

// ToDnaBases returns a slice of dna.Base that represents the same sequence
// of bases present in fragment
func ToDnaBases(fragment *ThreeBit) []dna.Base {
	return RangeToDnaBases(fragment, 0, fragment.Len)
}
