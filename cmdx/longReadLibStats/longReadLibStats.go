package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"sort"
)

func longReadLibStats(inFq string, readLengths string) {
	var lens []int
	var tot int
	var out *fileio.EasyWriter

	if readLengths != "" {
		out = fileio.EasyCreate(readLengths)
	}

	reads := fastq.GoReadToChan(inFq)

	for read := range reads {
		tot += len(read.Seq)
		lens = append(lens, len(read.Seq))
		if readLengths != "" {
			fileio.WriteToFileHandle(out, fmt.Sprintf("%d", len(read.Seq)))
		}
	}

	sort.Slice(lens, func(i, j int) bool {
		return lens[j] < lens[i]
	})

	n50, _ := fasta.CalculateN50L50(lens, tot/2)
	fmt.Printf("N50: %d", n50)

	if readLengths != "" {
		exception.PanicOnErr(out.Close())
	}
}

func main() {
	var readLengths *string = flag.String("readLengths", "", "provide a filename for all the read lengths to be written out to")
	flag.Parse()
	longReadLibStats(flag.Arg(0), *readLengths)
}
