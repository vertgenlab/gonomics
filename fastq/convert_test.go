package fastq

import (
	"testing"
)

var fq *Fastq = Read("testdata/10x.barcoded_test.fastq")[0]

func TestQualToString(t *testing.T) {
	if QualString(fq.Qual) != Uint8QualToString(fq.Qual) {
		t.Errorf("Error: conversion was not successful...\n")
	}
}

func BenchmarkStringBuilder(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		QualString(fq.Qual)
	}
}

func BenchmarkFmtPrint(b *testing.B) {
	b.ReportAllocs()
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Uint8QualToString(fq.Qual)
	}
}

// Moved to benchmark after much improved performance in QualString
func Uint8QualToString(qual []uint8) string {
	var answer []rune = make([]rune, len(qual))
	for i := 0; i < len(qual); i++ {
		answer[i] = rune(qual[i] + asciiOffset)
	}
	return string(answer)
}
