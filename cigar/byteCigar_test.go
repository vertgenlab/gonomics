package cigar

import (
	"testing"
)

// Current declaration of the cigar struct
var c1 Cigar = Cigar{RunLength: 35, Op: 'M'}
var c2 Cigar = Cigar{RunLength: 2, Op: 'I'}
var c3 Cigar = Cigar{RunLength: 16, Op: 'D'}

// Byte Cigar light weight cigar
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
	threeFiveM := uint32(0) | uint32(35)<<4
	twoI := uint32(1) | uint32(2)<<4
	sixteenD := uint32(2) | uint32(16)<<4
	var bamCigar []uint32 = []uint32{threeFiveM, twoI, sixteenD}
	byteCigs := Uint32ToByteCigar(bamCigar)
	if string(ByteCigarToString(byteCigs)) != "35M2I16D" {
		t.Errorf("Error: could not convert %v to cigar...\n", byteCigs)
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
		ToString([]*Cigar{&c1, &c2, &c3})
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
