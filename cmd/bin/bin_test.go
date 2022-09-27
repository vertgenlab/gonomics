package main

import (
	"testing"
	"os"
	"log"
)

var (
	genome = "testdata/testGenome.fa"
)

func TestBinGenome(t *testing.T){
	binGenome(genome, 0, "testdata", 6)

	err1 := os.Remove("testdata/chr1.fa")
	err2 := os.Remove("testdata/_chr2_chr3.fa")
	err3 := os.Remove("testdata/chr4.fa")

	if err1 != nil || err2 != nil || err3 != nil {
		log.Panic(err1, err2, err3)
	}
}