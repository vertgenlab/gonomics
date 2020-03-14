package simpleGraph

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
	"testing"
)

func TestAlignSplit(t *testing.T) {
	var tileSize int = 32
	var numberOfReads int = 1000
	var readLength int = 150
	var mutations int = 0
	var numWorkers int = 8
	genome, _ := Read("testdata/gasAcu1.fa")

	simReads := RandomPairedReads(genome.Nodes, readLength, numberOfReads, mutations)
	fastq.WritePair("testdata/simReads_R1.fq", "testdata/simReads_R2.fq", simReads)
	log.Printf("Reading in the genome (simple graph)...\n")
	filenames := []string{"testdata/chrI.fa", "testdata/chrII.fa", "testdata/chrIII.fa", "testdata/chrIV.fa", "testdata/chrIX.fa", "testdata/chrM.fa", "testdata/chrUn.fa", "testdata/chrV.fa", "testdata/chrVI.fa", "testdata/chrVII.fa", "testdata/chrVIII.fa", "testdata/chrX.fa", "testdata/chrXI.fa", "testdata/chrXII.fa", "testdata/chrXIII.fa", "testdata/chrXIV.fa", "testdata/chrXIX.fa", "testdata/chrXV.fa", "testdata/chrXVI.fa", "testdata/chrXVII.fa", "testdata/chrXVIII.fa", "testdata/chrXX.fa", "testdata/chrXXI.fa"}
	//var chrGraphs []*SimpleGraph
	for i := 0; i < len(filenames); i++ {
		//chr, sizes := Read(filenames[i])
		//header := sam.ChromInfoMapSamHeader(sizes)

		tmpFile := fmt.Sprintf("tmp%d", i)
		GSWsBatchPair(filenames[i], "testdata/simReads_R1.fq", "testdata/simReads_R2.fq", tmpFile, numWorkers, tileSize)
		sorted, _ := sam.Read(tmpFile)
		sam.SortByName(sorted.Aln)
		samName := fmt.Sprintf("testdata/chr%d_Split.sam", i)
		sam.Write(samName, sorted)
		os.Remove(tmpFile)
	}
}
