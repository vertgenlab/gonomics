package main

import "testing"

var KmerSpectrumTests = []struct {
	FaFile       string
	BedFile      string
	OutFile      string
	ExpectedFile string
	KmerLength   int
}{
	{FaFile: "testdata/rand.fa",
		BedFile:      "testdata/regions.bed",
		OutFile:      "testdata/test.txt",
		ExpectedFile: "testdata/expected.txt",
		KmerLength:   3,
	},
}

func TestKmerSpectrum(t *testing.T) {
	var s Settings
	for _, v := range KmerSpectrumTests {
		s = Settings{
			FaFile:     v.FaFile,
			BedFile:    v.BedFile,
			OutFile:    v.OutFile,
			KmerLength: v.KmerLength,
		}
		kmerSpectrum(s)
	}
}
