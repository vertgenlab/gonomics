package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestSeek(t *testing.T) {
	sr := NewSeeker("testdata/cornerCases.fa", "")
	if dna.BasesToString(SeekByName(sr, "beta", 0, 2)) != "AG" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "beta", 1, 5)) != "GCCG" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "beta", 2, 16)) != "CCGA" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "gamma", 0, 100)) != "ACGT" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "gamma", 2, 2)) != "" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "gamma", 200, 2000)) != "" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "delta", 2, 3)) != "G" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "delta", 3, 4)) != "C" {
		t.Error("problem with fasta seek")
	}
	if dna.BasesToString(SeekByName(sr, "delta", 0, 160)) != "ACGCCCTAGT" {
		t.Error("problem with fasta seek")
	}
}
