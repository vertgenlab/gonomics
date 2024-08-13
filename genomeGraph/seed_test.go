package genomeGraph

import (
	"testing"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
)

// Benchmark Function
func BenchmarkExtendToTheRightDev(b *testing.B) {
    // Setup (Create representative Node and FastqBig instances)
	genome := Read("testdata/mini.gg")
    fq := fastq.ToFastqBig(fastq.Fastq{Name: "read", Seq: dna.StringToBases("CGTCATGTGCTACTGAATGCTGTACTCGTTACACGTTGT")})
    // Benchmark Loop
    for i := 0; i < b.N; i++ {
        _ = extendToTheRightDev(&genome.Nodes[3], &fq, 0, 0, true, nil)
    }
	b.ReportAllocs()
}