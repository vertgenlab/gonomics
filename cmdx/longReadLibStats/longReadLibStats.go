package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sort"
)

func longReadLibStats(inFq string, outfile string, readLengths string) {
	var lens []int
	var tot int
	var out *fileio.EasyWriter
	var err error

	if readLengths != "" {
		out = fileio.EasyCreate(readLengths)
	}

	reads := fastq.GoReadToChan(inFq)
	outStats := fileio.EasyCreate(outfile)

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

	fileio.WriteToFileHandle(outStats, fmt.Sprintf("Total number of reads: %d", len(lens)))
	fileio.WriteToFileHandle(outStats, fmt.Sprintf("N50: %d", n50))

	if readLengths != "" {
		err = out.Close()
		exception.PanicOnErr(err)
	}

	err = outStats.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("longReadLibStats -- Provide statistics for a long read fastq dataset including N50 and number of reads\n" +
		"N50 is length of the read at which half the bases in the library exist at read lengths larger than the N50 read\n" +
		"longReadLibStats [options] in.fq out.stats.txt\n" +
		"Options:\n")
	flag.PrintDefaults()
}

func main() {
	var readLengths *string = flag.String("readLengths", "", "provide a filename for all the read lengths to be written out to")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != 2 {
		flag.Usage()
		log.Fatalf("Expecting 2 agument but got %d", len(flag.Args()))
	}

	longReadLibStats(flag.Arg(0), flag.Arg(1), *readLengths)
}
