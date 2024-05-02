package cigar

import (
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
	perfect := []ByteCigar{
		{RunLen: 150, Op: 'M'},
	}
	result0 := SoftClipBases(0, 150, perfect)

	expected0 := []ByteCigar{
		{RunLen: 150, Op: 'M'},
	}

	if ByteCigarToString(expected0) != ByteCigarToString(result0) {
		t.Errorf("Test case 1 failed. Expected %s, but got %s", ByteCigarToString(expected0), ByteCigarToString(result0))
	}

	cig := []ByteCigar{
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'D'},
		{RunLen: 5, Op: 'M'},
	}

	// Test case 1: front = 2, lengthOfRead = 20
	expected1 := []ByteCigar{
		{RunLen: 2, Op: 'S'},
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'D'},
		{RunLen: 5, Op: 'M'},
		{RunLen: 3, Op: 'S'},
	}
	result1 := SoftClipBases(2, 20, cig)

	if ByteCigarToString(result1) != ByteCigarToString(expected1) {
		t.Errorf("Test case 2 failed. Expected %s, but got %s", ByteCigarToString(expected1), ByteCigarToString(result1))
	}

	insertions := []ByteCigar{
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'I'},
		{RunLen: 5, Op: 'M'},
	}

	expected2 := []ByteCigar{
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'I'},
		{RunLen: 5, Op: 'M'},
		{RunLen: 5, Op: 'S'},
	}
	result2 := SoftClipBases(0, 25, insertions)
	if ByteCigarToString(result2) != ByteCigarToString(expected2) {
		t.Errorf("Test case 3 failed. Expected %s, but got %s", ByteCigarToString(expected2), ByteCigarToString(result2))
	}

	// Test case 3: front = 5, lengthOfRead = 20
	five := []ByteCigar{
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'D'},
		{RunLen: 1, Op: 'I'},
		{RunLen: 4, Op: 'M'},
	}

	expected3 := []ByteCigar{
		{RunLen: 5, Op: 'S'},
		{RunLen: 10, Op: 'M'},
		{RunLen: 5, Op: 'D'},
		{RunLen: 1, Op: 'I'},
		{RunLen: 4, Op: 'M'},
	}
	result3 := SoftClipBases(5, 20, five)

	if ByteCigarToString(result3) != ByteCigarToString(expected3) {
		t.Errorf("Test case 4 failed. Expected %s, but got %s", ByteCigarToString(expected3), ByteCigarToString(result3))
	}
}

// TestReverseByteCigars tests the reversing functionality of ReverseByteCigars function.
func TestReverseByteCigars(t *testing.T) {
	// Define a test case with an initial slice of ByteCigars.
	cigars := []ByteCigar{
		{RunLen: 1, Op: Match},
		{RunLen: 2, Op: Insertion},
		{RunLen: 3, Op: Deletion},
	}

	// What we expect after the reverse operation.
	expected := []ByteCigar{
		{RunLen: 3, Op: Deletion},
		{RunLen: 2, Op: Insertion},
		{RunLen: 1, Op: Match},
	}

	// Perform the reverse operation using the function under test.
	reversedCigars := ReverseBytesCigar(cigars)

	// Check if the reversed slice matches the expected slice.
	if ByteCigarToString(reversedCigars) != ByteCigarToString(expected) {
		t.Errorf("Error: ReverseByteCigars failed - %v != %v", reversedCigars, expected)
	}
}

func TestUint32Cigar(t *testing.T) {
	cigars := []CigarOp{36, 160, 34}
	expected := "2S10M2D"
	if Uint32ToString(cigars) != expected {
		t.Errorf("Error: Uint32Cigar conversion failed %s != %s...", Uint32ToString(cigars), expected)
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

func BenchmarkReverseCigar(b *testing.B) {
	cigars := []ByteCigar{
		{RunLen: 1, Op: Match},
		{RunLen: 2, Op: Insertion},
		{RunLen: 3, Op: Deletion},
		{RunLen: 1, Op: Match},
	}
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		ReverseBytesCigar(cigars)
	}
}
