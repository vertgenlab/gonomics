package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/sam"
	"testing"
)

func TestTagToReadGroup(t *testing.T) {
	in := "testdata/tagTest.bam"
	valueFile := "testdata/tagTestValues.txt"
	targetTag := "CB"
	out := "testdata/testOutput.bam"

	bamTagToReadGroup(in, out, targetTag, valueFile)

	inOutput, _ := sam.Read(in)
	testOutput, testHeader := sam.Read(out)
	truthOutput, truthHeader := sam.Read("testdata/output.bam")

	if len(testHeader.Text) != len(truthHeader.Text) {
		t.Error("problem making header")
	}

	for i := range testHeader.Text {
		if testHeader.Text[i] != truthHeader.Text[i] {
			t.Error("problem making header")
		}
	}

	if len(truthOutput) != len(inOutput) {
		t.Error("PROBLEM WITH OUTPUT TESTFILE, CHECK CAREFULLY")
	}

	if len(testOutput) != len(truthOutput) {
		t.Error("problem writing records")
	}

	var err error
	for i := range testOutput {
		err = sam.ParseExtra(&testOutput[i])
		exception.PanicOnErr(err)
		err = sam.ParseExtra(&truthOutput[i])
		exception.PanicOnErr(err)

		if testOutput[i].Extra != truthOutput[i].Extra {
			t.Error("problem writing extra")
		}
	}
}
