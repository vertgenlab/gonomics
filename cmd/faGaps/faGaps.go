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

func usage() {
	fmt.Print(
		"faGap - a program to investigate regions containing gaps\n" +
			"Usage:\n" +
			"  faGap [options] input.fa\n\n" +
			"Options:\n" +
			"  -b, --bed\toutput bed file containing regions outside gaps\n" +
			"  -s, --split\tbreak fasta into multiple records, split by gapped regions\n" +
			"  -o, --out\twrite to filename (default: /dev/stdout)\n\n")
}

var bedOut = flag.Bool("`'`b", false, "`--bed`output bed file containing regions outside gaps")
var faOut = flag.Bool("s", false, "`--split`break fasta into multiple records, split by gapped regions")
var outFile = flag.String("o", "/dev/stdout", "`--out`write to filename (default: /dev/stdout)")

func init() {
	flag.BoolVar(bedOut, "bed", false, "output bed file containing regions outside gaps")
	flag.BoolVar(faOut, "split", false, "break fasta into multiple records, split by gapped regions")
	flag.StringVar(outFile, "out", "/dev/stdout", "write to filename (default: /dev/stdout)")
	flag.Usage = usage
}

func main() {
	var expectedNumArgs int = 1
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	var inFile string = flag.Arg(0)
	if *bedOut {
		faNoGap(inFile, *outFile)
	} else if *faOut {
		faSplitByNs(inFile, *outFile)
	} else {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
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

func faSplitByNs(filename string, outFile string) {
	fa := fasta.Read(filename)
	out := fileio.EasyCreate(outFile)
	for i := 0; i < len(fa); i++ {
		fasta.WriteToFileHandle(out, chrSplitByNs(fa[i]), 50)
	}
}
