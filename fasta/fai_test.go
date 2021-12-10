package fasta

import (
	"testing"
)

func TestCreateIndex(t *testing.T) {
	actual := CreateIndex("testdata/testOne_fmt.fa")
	expected := readIndex("testdata/testOne_fmt.fa.fai")
	if actual.String() != expected.String() {
		t.Errorf("gonomics and samtools fasta indexes do not match.\ngonomics:\n%s\n\nsamtools:\n%s\n", actual, expected)
	}

	actual = CreateIndex("testdata/cornerCases.fa")
	expected = readIndex("testdata/cornerCases.fa.fai")
	if actual.String() != expected.String() {
		t.Errorf("gonomics and samtools fasta indexes do not match.\ngonomics:\n%s\n\nsamtools:\n%s\n", actual, expected)
	}
}
