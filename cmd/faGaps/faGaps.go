package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func faNoGap(inFile string, outFile string) {
	records := fasta.Read(inFile)
	beds := bed.UngappedRegionsAllFromFa(records)
	bed.Write(outFile, beds, 3)
}

func faSplitByNs(filename string, outFile string) {
	reader := fileio.EasyOpen(filename)
	writer := fileio.EasyCreate(outFile)
	defer reader.Close()
	defer writer.Close()
	for fa, done := fasta.NextFasta(reader); !done; fa, done = fasta.NextFasta(reader) {
		fasta.WriteToFileHandle(writer, chrSplitByNs(fa), 50)
	}
}

func usage() {
	fmt.Print(
		"faGap - a program to investigate regions containing gaps\n\n" +
			"Usage:\n" +
			"  faGap [options] in.fa out.file\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	var expectedNumArgs int = 2
	var bedOut = flag.Bool("bed", false, "find genomic coordinates containing regions outside gaps `.bed`")
	var faOut = flag.Bool("fasta", false, "split fasta into several records using gapped coordinates `.fa/.gz`\n")

	flag.Parse()

	var inFile string = flag.Arg(0)
	var outFile string = flag.Arg(1)

	if *bedOut {
		faNoGap(inFile, outFile)
	} else if *faOut {
		faSplitByNs(inFile, outFile)
	} else {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	}
}

//TODO: Belongs in the fasta package, but bed and fasta becomes an illegal import cycle
func chrSplitByNs(chr *fasta.Fasta) []*fasta.Fasta {
	unGapped := bed.UngappedRegionsFromFa(chr)
	var answer []*fasta.Fasta = make([]*fasta.Fasta, len(unGapped))
	for i := 0; i < len(unGapped); i++ {
		answer[i] = &fasta.Fasta{Name: fmt.Sprintf("%s_%d_%d", unGapped[i].Chrom, unGapped[i].ChromStart, unGapped[i].ChromEnd), Seq: chr.Seq[unGapped[i].ChromStart:unGapped[i].ChromEnd]}
	}
	return answer
}
