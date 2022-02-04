// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfToFa(inVcfFilename string, inFaFilename string, outFaFilename string, useAlt bool) {
	var seqsOrdered []fasta.Fasta
	var seqsLookup fasta.FastaMap
	var vcfRecords <-chan vcf.Vcf

	// We hang on to the ordered version of the fasta sequences so that we could
	// output the edited sequences (those in the map) in the same order that
	// we originally read them, but we do not have that currently
	seqsOrdered = fasta.Read(inFaFilename)
	seqsLookup = fasta.ToMap(seqsOrdered)

	// we could do some sort of header check here to make sure the fa matches the vcf
	vcfRecords, _ = vcf.GoReadToChan(inVcfFilename)
	for v := range vcfRecords {
		if !(vcf.IsBiallelic(v) && vcf.IsSubstitution(v)) {
			log.Fatal("Error: currently we only handle biallelic substitutions\n")
		}

		// check to make sure the fa matches the vcf before we edit
		if seqsLookup[v.Chr][v.Pos-1] != dna.StringToBase(v.Ref) {
			log.Fatal("Error: base in fasta didn't match ref base from VCF record\n")
		}

		if useAlt {
			seqsLookup[v.Chr][v.Pos-1] = dna.StringToBase(v.Alt[0])
		}
	}

	fasta.WriteMap(outFaFilename, seqsLookup)
}

func usage() {
	fmt.Print(
		"vcfToFa - Use the variant data in the vcf to edit an input fasta of the reference.\n\n" +
			"Usage:\n" +
			"  vcfToFa [options] input.vcf input.fa output.fa\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var useAlt *bool = flag.Bool("useAlt", false, "Retains only biallelic variants in the output file.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inVcfFilename := flag.Arg(0)
	inFaFilename := flag.Arg(1)
	outFaFilename := flag.Arg(2)

	vcfToFa(inVcfFilename, inFaFilename, outFaFilename, *useAlt)
}
