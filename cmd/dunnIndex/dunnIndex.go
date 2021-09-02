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

func dunnIndex(bedFile string, alnFile string, groupFileName string, outFile string) {
	b := bed.Read(bedFile)
	aln := fasta.Read(alnFile)
	g := popgen.ReadGroups(groupFileName)

	for i := 0; i < len(b); i++ {
		b[i].Annotation = make([]string, 3)
		dunn, S, missing := popgen.Dunn(b[i], aln, g)
		b[i].Annotation[0] = fmt.Sprintf("%f", dunn)
		b[i].Annotation[1] = fmt.Sprintf("%v", S)
		b[i].Annotation[2] = missing
	}
	bed.Write(outFile, b)
}

func usage() {
	fmt.Print(
		"dunnIndex - Computes the Dunn Index based on variable SNPs for each input bed region of a multiple alignment.\n" +
			"Groups should be specified in a .list file with group names on lines beginning with a 'greater than' symbol and all group members on the following lines.\n" +
			"Returns a bed file with the Dunn Index, number of segregating sites, and missing individuals in the annotation columns as tab separated fields.\n" +
			"Usage:\n" +
			"dunnIndex regions.bed aln.multi.fa  groups.list outfile.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
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
	groupFile := flag.Arg(2)
	outFile := flag.Arg(3)

	dunnIndex(bedFile, alnFile, groupFile, outFile)
}
