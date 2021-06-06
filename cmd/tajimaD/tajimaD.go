// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/popgen"
	"log"
)

func tajimaD(inFile string, alnFile string, groupFile string, outFile string) {
	b := bed.Read(inFile)
	aln := fasta.Read(alnFile)

	if groupFile == "" {
		for i := 0; i < len(b); i++ {
			b[i].Annotation = make([]string, 1)
			b[i].Annotation[0] = fmt.Sprintf("%f", popgen.TajimaFromBedNoGroup(b[i], aln))
		}
	} else {
		for i := 0; i < len(b); i++ {
			b[i].Annotation = make([]string, 2)
			tajima, missing := popgen.TajimaFromBed(b[i], aln, popgen.ReadGroups(groupFile))
			b[i].Annotation[0] = fmt.Sprintf("%f", tajima)
			b[i].Annotation[1] = missing
		}
	}
	bed.Write(outFile, b, 7)
}

func usage() {
	fmt.Print(
		"tajimaD - Computes Tajima's D for each region in an input bed for a given multiFa alignment.\n" +
			"Usage:\n" +
			"tajimaD regions.bed aln.multi.fa outfile.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var groupFile *string = flag.String("groups", "", "Restricts analysis to a subset of the alignment.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	bedFile := flag.Arg(0)
	alnFile := flag.Arg(1)
	outFile := flag.Arg(2)

	tajimaD(bedFile, alnFile, *groupFile, outFile)
}
