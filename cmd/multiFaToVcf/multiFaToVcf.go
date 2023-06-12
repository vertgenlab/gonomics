// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

func multiFaToVcf(inFile string, chr string, outFile string, substitutionsOnly bool, retainN bool) {
	f := fasta.Read(inFile)
	out := fileio.EasyCreate(outFile)
	header := vcf.NewHeader("")
	vcf.NewWriteHeader(out, header)
	convert.PairwiseFaToVcf(f, chr, out, substitutionsOnly, retainN)
	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"multiFaToVcf - Generates a VCF file from an input pairwise multiFa alignment with the first entry as the reference." +
			"Note that deletions in the first position of an alignment will not appear in the output Vcf.\n" +
			"Usage:\n" +
			"multiFaToVcf multi.Fa chromName out.vcf \n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var substitutionsOnly *bool = flag.Bool("substitutionsOnly", false, "Retain only substitutions in the output VCF.")
	var retainN *bool = flag.Bool("retainN", false, "By default, multiFaToVcf does not retain positions with N in the \n"+
		"alt or ref in the output VCF. This option overrides this behavior.")
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

	multiFaToVcf(inFile, chromName, outFile, *substitutionsOnly, *retainN)
}
