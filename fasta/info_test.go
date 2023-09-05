package fasta

import (
	"testing"

	"github.com/vertgenlab/gonomics/chromInfo"
)

func TestIsFasta(t *testing.T) {
	if !IsFasta("test.fa") {
		t.Errorf("Problem is IsFasta")
	}
	if !IsFasta("test.fa.gz") {
		t.Errorf("Problem is IsFasta")
	}
	if !IsFasta("test.fasta") {
		t.Errorf("Problem is IsFasta")
	}
	if !IsFasta("test.fasta.gz") {
		t.Errorf("Problem is IsFasta")
	}
	if IsFasta("test.vcf") {
		t.Errorf("Problem is IsFasta")
	}
	if IsFasta("test.vcf.gz") {
		t.Errorf("Problem is IsFasta")
	}
}

func TestToChromInfo(t *testing.T) {
	input := []Fasta{{"apple", seqOneA}, {"banana", seqOneB}, {"carrot", seqOneC}}
	output := ToChromInfo(input)

	expected := []chromInfo.ChromInfo{
		{Name: "apple", Size: len(seqOneA), Order: 0},
		{Name: "banana", Size: len(seqOneB), Order: 1},
		{Name: "carrot", Size: len(seqOneC), Order: 2},
	}

	if len(output) != len(expected) {
		t.Error("problem with ToChromInfo")
	}

	for i := range output {
		if output[i] != expected[i] {
			t.Error("problem with ToChromInfo")
		}
	}
}
