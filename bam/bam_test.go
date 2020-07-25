package bam

import (
	"github.com/vertgenlab/gonomics/sam"
	"testing"
)

//readBamTests are the files pairs used to test our binary reader functions
var readBamTests = []struct {
	bam string
	sam string
}{
	{"testdata/small.bam", "testdata/small.sam"},
	{"testdata/tenXbarcodeTest.bam", "testdata/tenXbarcodeTest.sam"},
}

//TestBamToSamReader will convert a bam file into a sam record and perform a comparison with the same sam file that was created using samtools view.
func TestBamToSamReader(t *testing.T) {
	for _, test := range readBamTests {
		bamFile := Read(test.bam)
		samFile, err := sam.Read(test.sam)
		if err != nil {
			t.Errorf("Error: There was a problem reading in the sam file...\n")
		}
		if len(bamFile) != len(samFile.Aln) {
			t.Errorf("Error: File lines are not equal...\n")
		}
		for i := 0; i < len(bamFile); i++ {
			if !sam.IsEqualDebug(bamFile[i], samFile.Aln[i]) {
				t.Fatalf("Error: Did not create the same sam file as samtools view...\n")
			}
		}
	}
}

//BenchmarkSamReader will benchmark the reading speed of a basic sam text file. (This is used to compare with reading a bam file)
func BenchmarkSamReader(b *testing.B) {
	var samFile *sam.Sam
	var err error
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		for _, test := range readBamTests {
			samFile, err = sam.Read(test.sam)
		}
	}
	if err != nil {
		b.Errorf("Error: There was a problem reading in the sam file after %d lines...\n", len(samFile.Aln))
	}
}

//BenchmarkBamReader will benchmark the speed of decoding a bam file and convert the data into sam records. (This is used to compare with reading a sam file)
func BenchmarkBamReader(b *testing.B) {
	var bamFile []*sam.SamAln
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		for _, test := range readBamTests {
			bamFile = Read(test.bam)
		}

	}
	if bamFile == nil {
		b.Errorf("Error: There was a problem reading in the sam file after %d lines...\n", len(bamFile))
	}
}
