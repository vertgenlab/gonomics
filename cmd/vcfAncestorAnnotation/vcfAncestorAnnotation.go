package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfAncestorAnnotation(inFile string, faFile string, outFile string) {
	ch, header := vcf.GoReadToChan(inFile)
	records := fasta.Read(faFile)
	out := fileio.EasyCreate(outFile)
	defer out.Close()

	vcf.NewWriteHeader(out, header)

	for v := range ch {
		vcf.VcfAnnotateAncestorFromFa(v, records)
		vcf.WriteVcf(out, v)
	}
}

func usage() {
	fmt.Print(
		"vcfAncestorAnnotation - Adds ancestral allele to the INFO column of entries in a VCF file.\n" +
			"The VCF file must be split by chromosome beforehand, as the fasta file is formatted with the first entry as the reference sequence\n" +
			"and the second entry as the inferred ancestor sequence for a particular chromosome.\n" +
			"Usage:\n" +
			"vcfAncestorAnnotation input.vcf alignment.fa output.vcf\n" +
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
	faFile := flag.Arg(1)
	outFile := flag.Arg(2)

	vcfAncestorAnnotation(inFile, faFile, outFile)
}
