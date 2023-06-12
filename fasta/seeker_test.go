package fasta

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
)

func TestSeek(t *testing.T) {
	sr := NewSeeker("testdata/cornerCases.fa", "")
	ans, err := SeekByName(sr, "beta", 0, 2)
	if err != nil {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "AG" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "beta", 1, 5)
	if err != nil {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "GCCG" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "beta", 2, 16)
	if err != nil && err != ErrSeekEndOutsideChr {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "CCGA" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "gamma", 0, 100)
	if err != nil && err != ErrSeekEndOutsideChr {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "ACGT" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "gamma", 2, 2)
	if err != nil {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "gamma", 200, 2000)
	if err != nil && err != ErrSeekStartOutsideChr {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "delta", 2, 3)
	if err != nil {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "G" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "delta", 3, 4)
	if err != nil {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "C" {
		t.Error("problem with fasta seek")
	}

	ans, err = SeekByName(sr, "delta", 0, 160)
	if err != nil && err != ErrSeekEndOutsideChr {
		t.Error(err)
	}
	if dna.BasesToString(ans) != "ACGCCCTAGT" {
		t.Error("problem with fasta seek")
	}
}
