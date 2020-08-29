package simpleGraph

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	//"sort"
)

type GirafGsw struct {
	ReadOne giraf.Giraf
	ReadTwo giraf.Giraf
}

func readFastqGsw(fileOne string, fileTwo string, answer chan<- fastq.PairedEndBig) {
	readOne, readTwo := fileio.NewSimpleReader(fileOne), fileio.NewSimpleReader(fileTwo)
	for fq, done := fastq.ReadFqBigPair(readOne, readTwo); !done; fq, done = fastq.ReadFqBigPair(readOne, readTwo) {
		answer <- *fq
	}
	close(answer)
}
