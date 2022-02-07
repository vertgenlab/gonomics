package dna

import (
	"log"
	"strings"
	"testing"
	"unicode/utf8"
)

var equalStringsAndBases = []struct {
	characters string // first input
	bases      []Base // second input
}{
	{"ACGTNacgtn", []Base{A, C, G, T, N, LowerA, LowerC, LowerG, LowerT, LowerN}},
	{"ACACGGGGNNNNNN", []Base{A, C, A, C, G, G, G, G, N, N, N, N, N, N}},
	{"NNNNnnnnTGTGTG", []Base{N, N, N, N, LowerN, LowerN, LowerN, LowerN, T, G, T, G, T, G}},
	{"NNNNgtgtacaNNN", []Base{N, N, N, N, LowerG, LowerT, LowerG, LowerT, LowerA, LowerC, LowerA, N, N, N}},
	{"nnnacgttgcannnaacc", []Base{LowerN, LowerN, LowerN, LowerA, LowerC, LowerG, LowerT, LowerT, LowerG, LowerC, LowerA, LowerN, LowerN, LowerN, LowerA, LowerA, LowerC, LowerC}},
}

func TestStringToBases(t *testing.T) {
	for _, test := range equalStringsAndBases {
		actual := StringToBases(test.characters)
		if CompareSeqsCaseSensitive(actual, test.bases) != 0 {
			t.Errorf("StringToBases(%s): expected %v, actual %v", test.characters, test.bases, actual)
		}
	}
}

func TestBasesToString(t *testing.T) {
	for _, test := range equalStringsAndBases {
		actual := BasesToString(test.bases)
		if strings.Compare(actual, test.characters) != 0 {
			t.Errorf("BasesToString(%v): expected %s, actual %v", test.bases, test.characters, actual)
		}
	}
}

var equalStringAndBase = []struct {
	character string
	base      Base
}{
	{"A", A},
	{"C", C},
	{"G", G},
	{"T", T},
	{"N", N},
	{"a", LowerA},
	{"c", LowerC},
	{"g", LowerG},
	{"t", LowerT},
	{"n", LowerN},
}

func TestStringToBase(t *testing.T) {
	for _, test := range equalStringAndBase {
		actual := StringToBase(test.character)
		if actual != test.base {
			t.Errorf("StringToBase(%s): expected %d, actual %d", test.character, test.base, actual)
		}
	}
}

func TestBaseToString(t *testing.T) {
	for _, test := range equalStringAndBase {
		actual := BaseToString(test.base)
		if actual != test.character {
			t.Errorf("BaseToString(%d): expected %s, actual %s", test.base, test.character, actual)
		}
	}
}

// BENCHMARK RESULTS
// BenchmarkBasesToStringViaSwitch-4         313066              3900 ns/op
// BenchmarkBasesToStringViaArray-4         1543102               763 ns/op
// BenchmarkStringToBasesViaSwitch-4         197355              5905 ns/op
// BenchmarkStringToBasesViaMap-4             32530             36990 ns/op
// BenchmarkStringToBasesViaArray-4          621732              1958 ns/op

// NOTE: converting from string to bases is significantly faster with an array, but i feel
// the implementation is quite messy, with little room for efficient error checking
// so I am leaving the current switch function, but people should be aware that for
// highly performance sensitive tasks, an array may be faster than the ByteToBase function.

func BenchmarkBasesToStringViaSwitch(b *testing.B) {
	seq := StringToBases(gdf2mrna)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for j := range seq {
			_ = BaseToRune(seq[j])
		}
	}
}

func BenchmarkBasesToStringViaArray(b *testing.B) {
	seq := StringToBases(gdf2mrna)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for j := range seq {
			_ = baseToStringArray[seq[j]]
		}
	}
}

func BenchmarkStringToBasesViaSwitch(b *testing.B) {
	for i := 0; i < b.N; i++ {
		stringToBasesSwitch(gdf2mrna)
	}
}

func BenchmarkStringToBasesViaMap(b *testing.B) {
	for i := 0; i < b.N; i++ {
		stringToBasesMap(gdf2mrna)
	}
}

func BenchmarkStringToBasesViaArray(b *testing.B) {
	arr := makeStringToBaseArray()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		stringToBasesArray(gdf2mrna, arr)
	}
}

func stringToBasesSwitch(s string) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index := range s {
		answer[index] = byteToBase(s[index])
	}
	return answer
}

func stringToBasesMap(s string) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index := range s {
		answer[index] = stringBaseMap[s[index]]
	}
	return answer
}

func stringToBasesArray(s string, arr []Base) []Base {
	answer := make([]Base, utf8.RuneCountInString(s))

	for index := range s {
		answer[index] = arr[s[index]]
	}
	return answer
}

var baseToStringArray = []byte{'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', '-', '.', '*'}

func makeStringToBaseArray() []Base {
	var stringToBaseArray = make([]Base, 't'+1)
	stringToBaseArray['A'] = A
	stringToBaseArray['C'] = C
	stringToBaseArray['G'] = G
	stringToBaseArray['T'] = T
	stringToBaseArray['N'] = N
	stringToBaseArray['a'] = LowerA
	stringToBaseArray['c'] = LowerC
	stringToBaseArray['g'] = LowerG
	stringToBaseArray['t'] = LowerT
	stringToBaseArray['n'] = LowerN
	stringToBaseArray['-'] = Gap
	stringToBaseArray['.'] = Dot
	stringToBaseArray['*'] = Nil
	return stringToBaseArray
}

var stringBaseMap = map[byte]Base{
	'A': A,
	'C': C,
	'G': G,
	'T': T,
	'a': LowerA,
	'c': LowerC,
	'g': LowerG,
	't': LowerT,
	'N': N,
	'n': LowerN,
	'-': Gap,
	'.': Dot,
	'*': Nil,
}

func byteToBase(b byte) Base {
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
	case '*':
		return Nil
	case '.':
		return Dot
	default:
		log.Fatalf("Error: unexpected character in dna %c\n", b)
		return N
	}
}
