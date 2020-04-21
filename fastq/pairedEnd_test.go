package fastq

import (
	"fmt"
	"testing"
)

func TestReadPairChan(t *testing.T) {
	fastqPipe := make(chan *PairedEnd, 824)
	var readOne string = "testdata/simReads_R1.fq"
	var readTwo string = "testdata/simReads_R2.fq"
	go PairEndToChan(readOne, readTwo, fastqPipe)
	alpha := ReadPairs(readOne, readTwo)
	beta := make([]*PairedEnd, 0)
	for readPair := range fastqPipe {
		beta = append(beta, readPair)
	}
	if len(alpha) != len(beta) {
		fmt.Errorf("Error: Read fastq pair to channel did not result in the same length...\n")
	}
}
