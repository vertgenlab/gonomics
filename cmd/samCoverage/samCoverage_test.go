package main

import (
	"testing"
)

func TestSamCoverage(t *testing.T) {
	samCoverage("testdata/RandReads.sorted.sam", "testdata/test.txt", false)
}
