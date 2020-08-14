package main 

import(
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
)

func multiFaToVcf(inFile string, chr string, outFile string) {
	f := fasta.Read(inFile)
	v := convert.PairwiseFaToVcf(f, chr)
	vcf.Write(outFile, v)
}

func usage() {
	fmt.Print(
		"multiFaToVcf - Generates a VCF file from an input pairwise multiFa alignment with the first entry as the reference.\n" +
		"Usage:\n" +
		"multiFaToVcf multi.Fa chromName out.vcf \n" +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	inFile := flag.Arg(0)
	chromName := flag.Arg(1)
	outFile := flag.Arg(2)

	multiFaToVcf(inFile, chromName, outFile)
}