package fasta

import "testing"

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
