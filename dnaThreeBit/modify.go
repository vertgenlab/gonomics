package dnaThreeBit

import (
	"bytes"
	"log"
	"unicode/utf8"
)

// Append adds "b" to the end of "fragment."
// "fragment" can be nil.
func Append(fragment *ThreeBit, b *ThreeBitBase) *ThreeBit { // too confusing to have Append and append in this package?
	if fragment == nil {
		return &ThreeBit{Seq: uint64{b << 61}, Len: 1}
	}
	var basesInLastIdx int = fragment.Len % 21
	if basesInLastIdx == 0 { // there is no more room in the slice of uint64
		fragment.Seq = append(fragment.Seq, b << 61) // shift the new base all the way to the left
	} else { // there is room for at least one more base in the last uint64
		var lastIdx int = len(fragment.Seq) - 1
		fragment.Seq[lastIdx] = fragment.Seq[lastIdx] | (b << (61 - basesInLastIdx * 3)) // bit-wise or adds the new base
	}
	fragment.Len++
	return fragment
}

// Cat appends "b" to "a."  "a" is changed to be both of them, and "b" is unchanged.
// It is quickest to have "a" be the longer sequence.
func Cat(a *ThreeBit, b *ThreeBit) { // I was tempted to call this ligate
	if b == nil {
		return
	}
	for i := 0; i < b.Len; i++ {
		a = Append(a, getThreeBitBase(b, i))
	}
}

func Copy(a *ThreeBit) *ThreeBit {
	if a == nil {
		return nil
	}
	answer := &ThreeBit{Seq: make(uint64[], a.Len), Len: a.Len}
	copy(answer.Seq, a.Seq)
	return answer
}

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

func FromString(s string) *ThreeBit {
	answer := &ThreeBit{Seq:[]uint64{}, Len:0}
	for _, runeValue := range s {
		answer = Append(answer, RuneToThreeBitBase(runeValue))
	}
	return answer
}

func AppendBytes(a *ThreeBit, b []byte) *ThreeBit {
        if a == nil {
                a = &ThreeBit{
        }
        var basesInLastIdx int = fragment.Len % 21
        if basesInLastIdx == 0 { // there is no more room in the slice of uint64
                fragment.Seq = append(fragment.Seq, b << 61) // shift the new base all the way to the left
        } else { // there is room for at least one more base in the last uint64
                var lastIdx int = len(fragment.Seq) - 1
                fragment.Seq[lastIdx] = fragment.Seq[lastIdx] | (b << (61 - basesInLastIdx * 3)) // bit-wise or adds the new base
        }
        fragment.Len++
        return fragment
}

func FromBytes(b []byte) *ThreeBit {


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

