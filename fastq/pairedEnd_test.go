package fastq

import (
	"testing"
)

//TODO: files to read was modified because others didn't exist
// would be good to make this a tougher test
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
		t.Errorf("Error: Read fastq pair to channel did not result in the same length...\n")
	}
}
