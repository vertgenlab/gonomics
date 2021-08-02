package main

import (
	"io/ioutil"
	"os"
	"testing"
)

// minimal test for now until more complex functionality is implemented
// correctness is tested in the Sam package.
func TestPileup(t *testing.T) {
	tmpfile, _ := ioutil.TempFile("", "")
	pileup("../../sam/testdata/peak.bam", tmpfile.Name(), Settings{minDp: 0})
	tmpfile.Close()
	os.Remove(tmpfile.Name())
}
