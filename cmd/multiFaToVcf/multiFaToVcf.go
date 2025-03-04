// Command Group: "Data Conversion"

// Generates a VCF file from an input pairwise or three-way multiFa alignment with the first entry as the reference
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

func multiFaToVcf(inFile string, chr string, outFile string, substitutionsOnly bool, retainN bool, secondQueryName string) {
	f := fasta.Read(inFile)
	out := fileio.EasyCreate(outFile)
	header := vcf.NewHeader()
	vcf.NewWriteHeader(out, header)

	if secondQueryName != "" { // if a secondQueryName is specified
		// first check that secondQueryName can be found in the multiFa
		var fMap = fasta.ToMap(f)
		var secondQuerySequence, found = fMap[secondQueryName]
		if !found {
			log.Fatalf("Error: second query name is specified, but not found in the input multiFa file.\n")
		}
		// get the secondQuery Fasta, i.e. getting both Name and Seq
		var secondQueryFasta = fasta.Fasta{Name: secondQueryName, Seq: secondQuerySequence}
		// make an updated, edited multiFa, with the 2 sequences to convert to vcf in pairwise mode
		var fEdited = []fasta.Fasta{f[0], secondQueryFasta}
		// call convert.PairwiseFaToVcf on the edited multiFa
		convert.PairwiseFaToVcf(fEdited, chr, out, substitutionsOnly, retainN)
	} else {
		if len(f) == 2 {
			convert.PairwiseFaToVcf(f, chr, out, substitutionsOnly, retainN)
		} else if len(f) == 3 {
			convert.ThreeWayFaToVcf(f, chr, out)
		} else {
			log.Fatalf("Error: expecting 2 or 3 sequences in the input FASTA.\n")
		}
	}

	err := out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"multiFaToVcf - Generates a VCF file from an input multiFa alignment with the first entry as the reference.\n" +
			"If the input multiFa is pairwise alignment, checks for substitutions as well as indels,\n" +
			"but if the input multiFa is three-way alignment, it only checks for substitutions.\n" +
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
	var secondQueryName *string = flag.String("secondQueryName", "", "Specify the name of the second query sequence. The first query sequence is assumed to be the reference sequence, i.e. the first sequence in the multiFa.")
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

	multiFaToVcf(inFile, chromName, outFile, *substitutionsOnly, *retainN, *secondQueryName)
}
