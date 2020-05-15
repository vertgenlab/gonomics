package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfFormat(infile string, outfile string, ensemblToUCSC bool, UCSCToEnsembl bool, fixVcfRecords bool, ref string) {
	ch := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	defer out.Close()

	if ensemblToUCSC && UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatable.")
	}

	if fixVcfRecords {
		ref := fasta.FastaToMap(fasta.Read(ref))
		if ensemblToUCSC {
			for v := range ch {
				vcf.FixVcf(v, ref)
				v.Chr = convert.EnsemblToUCSC(v.Chr)
				vcf.WriteVcf(out.File, v)
			}
		} else if UCSCToEnsembl {
			for v := range ch {
				vcf.FixVcf(v, ref)
				v.Chr = convert.UCSCToEnsembl(v.Chr)
				vcf.WriteVcf(out.File, v)
			}
		}
	} else {
		if ensemblToUCSC {
			for v := range ch {
				v.Chr = convert.EnsemblToUCSC(v.Chr)
				vcf.WriteVcf(out.File, v)
			}
		} else if UCSCToEnsembl {
			for v := range ch {
				v.Chr = convert.UCSCToEnsembl(v.Chr)
				vcf.WriteVcf(out.File, v)
			}
		}
	}
}

func usage() {
	fmt.Print(
		"vcfFormat: Options alter VCF formatting.\n" +
			"Usage:\n" +
			"vcfFormat input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}


func main() {
	var expectedNumArgs int = 2
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")
	var fixVcfRecords *bool = flag.Bool("fix", false, "Fixes improperly formatted vcf records (e.g. '-' in ALT field")
	var ref *string = flag.String("ref", "", "Reference fasta. Only needed if using -fix.")

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

	vcfFormat(infile, outfile, *ensemblToUCSC, *UCSCToEnsembl, *fixVcfRecords, *ref)
}
