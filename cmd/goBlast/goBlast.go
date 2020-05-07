package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"log"
	"os"
)

func usage() {
	fmt.Print(
		"goBlast - quick tools for using the ncbi database\n" +
			"Usage:\n" +
			"goBlast [options] args\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	if len(os.Args) < 2 {
		usage()
	}
	var split *bool = flag.Bool("split", false, "contigs.fa prefix number\ndivide fasta into n files")
	var search *bool = flag.Bool("search", false, "summary.txt taxonId.txt *.tsv\nhard filter records reads that match given taxIds")
	var filter *bool = flag.Bool("filter", false, "summary.txt merge.fasta.gz *.fa\nfilter fasta and merge")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.LstdFlags)
	flag.Parse()
	if *split {
		EfficientTranspose(flag.Arg(0), flag.Arg(1), common.StringToInt(flag.Arg(2)))
	} else if *search {
		cmdArg := flag.Args()
		WrapAllDbFiles(cmdArg[0], cmdArg[1], cmdArg[2:])
	} else if *filter {
		cmdArg := flag.Args()
		FinalContigsToCanu(cmdArg[0], cmdArg[1], cmdArg[2:])
	} else {
		usage()
		flag.PrintDefaults()
	}
}
/*
func gatherFastafiles(output string, files []string) {
	faChan := make(chan *fasta.Fasta, 824)
	go fasta.MultiFileChan(files, faChan)
	var wg sync.WaitGroup
	merged := fileio.EasyCreate(output)
	defer merged.Close()
	var threads int = 6
	wg.Add(threads)
	var reading sync.WaitGroup
	reading.Add(len(files))
	for _, contig := range files {
		go fasta.ReadToChan(contig, faChan, &reading, false)
	}
	for i := 0; i < threads; i++ {
		go fasta.WritingChannel(merged, faChan, &wg)
	}
	reading.Wait()
	close(faChan)
}*/
