package dnaThreeBit

import (
	"bytes"
	"log"
	"github.com/vertgenlab/gonomics/dna"
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

func ThreeBitBaseToString(b ThreeBitBase) string {
	return string(ThreeBitBaseToRune(b))
}

func FromString(s string) *ThreeBit {
        answer := &ThreeBit{Seq:[]uint64{}, Len:0}
        for _, runeValue := range s {
                answer = Append(answer, RuneToThreeBitBase(runeValue))
        }
        return answer
}

func ToString(fragment *ThreeBit) string {
        var buffer bytes.Buffer

        for i:=0; i<fragment.Len; i++ {
                buffer.WriteRune(dna.BaseToRune(GetBase(fragment, i)))
        }
        return buffer.String()
}

func SectionToDnaBases(fragment *ThreeBit, start int, end int) []dna.Base {
	if end >= fragment.Len || start >= end {
		log.Fatalf("Error: unable to extract bases from %d to %d from a sequence of length %d\n", start, end, fragment.Len)
	}
        answer := make([]dna.Base, 0, end-start)
        for i := start; i < end; i++ {
                answer = append(answer, GetBase(fragment, i))
        }
        return answer
}

func ToDnaBases(fragment *ThreeBit) []dna.Base {
	return SectionToDnaBases(fragment, 0, fragment.Len)
}

