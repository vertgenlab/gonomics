// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

func vcfFormat(infile string, outfile string, ensemblToUCSC bool, UCSCToEnsembl bool, fixVcfRecords bool, ref string, clearInfo bool, tableOutput bool) {
	if ensemblToUCSC && UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatible.")
	}

	var maxAlts int
	if tableOutput {
		maxAlts = getMaxAltCount(infile)
	}

	ch, header := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	var err error

	var refMap map[string][]dna.Base
	var infoOrder []vcf.InfoHeader
	var formatOrder []vcf.FormatHeader

	if tableOutput {
		infoOrder, formatOrder = writeTableHeader(out, header, maxAlts)
	} else { // normal vcf output
		vcf.NewWriteHeader(out, header)
	}

	if fixVcfRecords {
		refMap = fasta.ToMap(fasta.Read(ref))
	}

	s := new(strings.Builder)
	for v := range ch {
		if clearInfo {
			v.Info = "."
		}
		if fixVcfRecords {
			v = vcf.FixVcf(v, refMap)
		}
		if ensemblToUCSC {
			v.Chr = convert.EnsemblToUCSC(v.Chr)
		}
		if UCSCToEnsembl {
			v.Chr = convert.UCSCToEnsembl(v.Chr)
		}
		if tableOutput {
			writeAsTable(s, out, v, header, infoOrder, formatOrder, maxAlts)
		} else { // normal vcf output
			vcf.WriteVcf(out, v)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"vcfFormat - Options alter VCF formatting.\n" +
			"Usage:\n" +
			"vcfFormat input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")
	var clearInfo *bool = flag.Bool("clearInfo", false, "Removes the information in the INFO field and replaces it with a '.'")
	var fixVcfRecords *bool = flag.Bool("fix", false, "Fixes improperly formatted vcf records (e.g. '-' in ALT field")
	var ref *string = flag.String("ref", "", "Reference fasta. Only needed if using -fix.")
	var tableOutput *bool = flag.Bool("csv", false, "Output as CSV file for spreadsheet analysis. Requires well-formed header.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	vcfFormat(infile, outfile, *ensemblToUCSC, *UCSCToEnsembl, *fixVcfRecords, *ref, *clearInfo, *tableOutput)
}
