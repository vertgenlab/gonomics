package cigar

import (
	"strings"
	"testing"
)

// Current declaration of the cigar struct.
var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}

// Byte Cigar light weight cigar.
var b1 ByteCigar = ByteCigar{RunLen: 35, Op: 'M'}
var b2 ByteCigar = ByteCigar{RunLen: 2, Op: 'I'}
var b3 ByteCigar = ByteCigar{RunLen: 16, Op: 'D'}

func TestBytesToCigar(t *testing.T) {
	var cigarbytes = []byte("35M2I16D")
	bc := ReadToBytesCigar(cigarbytes)
	if string(ByteCigarToString(bc)) != "35M2I16D" {
		t.Errorf("Error: could not convert %v to cigar...\n", bc)
	}
}

func TestUint32ToCigar(t *testing.T) {
	var query []uint32 = []uint32{560, 33, 258}
	byteCigs := Uint32ToByteCigar(query)
	if string(ByteCigarToString(byteCigs)) != "35M2I16D" {
		t.Errorf("Error: 35M!=560, 2I!=33, 16D!=258...\n")
	}
}

func TestToUint32(t *testing.T) {
	var answer []uint32 = []uint32{560, 33, 258}
	byteCigs := []ByteCigar{b1, b2, b3}
	code := ByteCigarToUint32(byteCigs)
	if len(answer) == len(code) {
		if code[0] != answer[0] || code[1] != answer[1] || code[2] != answer[2] {
			t.Errorf("Error: %d!=35M || %d!=2I || %d!=16D\n", code[0], code[1], code[2])
		}
	}
}

func TestSoftClipBases(t *testing.T) {
	testCases := []struct {
		name         string
		front        int
		lengthOfRead int
		cig          []ByteCigar
		expected     []ByteCigar
	}{
		{
			name:         "No Soft Clipping",
			front:        0,
			lengthOfRead: 10,
			cig:          []ByteCigar{{RunLen: 10, Op: 'M'}},
			expected:     []ByteCigar{{RunLen: 10, Op: 'M'}}, // No change
		},
		{
			name:         "Soft Clip at Front",
			front:        3,
			lengthOfRead: 10,
			cig:          []ByteCigar{{RunLen: 7, Op: 'M'}},
			expected:     []ByteCigar{{RunLen: 3, Op: 'S'}, {RunLen: 7, Op: 'M'}},
		},
		{
			name:         "Soft Clip at End",
			front:        0,
			lengthOfRead: 10,
			cig:          []ByteCigar{{RunLen: 5, Op: 'M'}},
			expected:     []ByteCigar{{RunLen: 5, Op: 'M'}, {RunLen: 5, Op: 'S'}},
		},
		{
			name:         "Soft Clip at Both Ends",
			front:        2,
			lengthOfRead: 10,
			cig:          []ByteCigar{{RunLen: 4, Op: 'M'}},
			expected:     []ByteCigar{{RunLen: 2, Op: 'S'}, {RunLen: 4, Op: 'M'}, {RunLen: 4, Op: 'S'}},
		},
		{
			name:         "Multiple CIGAR Ops",
			front:        1,
			lengthOfRead: 15,
			cig: []ByteCigar{
				{RunLen: 3, Op: 'M'},
				{RunLen: 2, Op: 'I'},
				{RunLen: 5, Op: 'D'},
			},
			expected: []ByteCigar{
				{RunLen: 1, Op: 'S'},
				{RunLen: 3, Op: 'M'},
				{RunLen: 2, Op: 'I'},
				{RunLen: 5, Op: 'D'},
				{RunLen: 9, Op: 'S'}, // 15 - 1 - 3 - 2 - 5 = 4
			},
		},
		// Add more test cases as needed to cover edge cases and different combinations
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			result := ByteCigarToString(SoftClipBases(tc.front, tc.lengthOfRead, tc.cig))
			if strings.Compare(result, ByteCigarToString(tc.expected)) != 0 {
				t.Errorf("Test Case %s: Expected %s, but got %s\n", tc.name, ByteCigarToString(tc.expected), result)
			}
		})
	}
}

func BenchmarkCigarBytesToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ByteCigarToString([]ByteCigar{b1, b2, b3})
	}
}

func BenchmarkCigarToString(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ToString([]Cigar{c1, c2, c3})
	}
}

func BenchmarkBytesToCigar(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var cigarbytes = []byte("35M2I16D")
	for n := 0; n < b.N; n++ {
		ReadToBytesCigar(cigarbytes)
	}
}

func BenchmarkStringToCigar(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	var cigarstring string = "35M2I16D"
	for n := 0; n < b.N; n++ {
		FromString(cigarstring)
	}
}
