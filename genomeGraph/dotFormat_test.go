package genomeGraph

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"log"
	"testing"
)

//TODO: consider dot format to graph
//     --> node_1            --> node_4
//node_0          --> node_3             --> node_6
//     --> node_2           -->  node_5

func TestDotFormat(t *testing.T) {
	gg := DotToGraph("testdata/dotFormat/seed_test.dot")
	var hippo *fastq.Fastq = &fastq.Fastq{Name: "hippoOne", Seq: dna.StringToBases("CCCATTCTTCCCTCCACAAGATCCAGATCCCCGTGTTTGGACCCCATAAGAATGAGAAGAATCCTCTCCGGTTTTGTTATCTGTTCTGAATGTGAAACTGGGATATCCACGTTGATCTCGGCTAAAATATAACAAACCTAAACAATTTGA"), Qual: fastq.ToQualUint8([]rune("JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ"))}
	seed := &SeedDev{TargetId: 3, TargetStart: 0, QueryStart: 23, Length: 10}
	hits := extendSeedTogether(seed, gg, hippo)
	printSeedDev(hits)
	log.Printf("Number of expanded seeds is: %d", len(hits))

}
