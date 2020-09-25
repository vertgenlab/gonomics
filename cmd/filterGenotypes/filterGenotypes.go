package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"filterGenotypes - filter genotyped Vcfs containing at leaste 3 samples with SNP regions where parental genomes are homozygous (and different) and have a clear heterozygous F1 from the parents\n\n" +
			"Usage:\n" +
			"  ./filterGenotypes [options] input.vcf output.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	var f1Genome *string = flag.String("f1", "", "F1 hybrid sample that appears heterozygous in genotype Vcf``")
	var sampleName *bool = flag.Bool("samples", false, "Get names of samples that appear in Vcf header. (Default: /dev/stdout)")

	var parentOne *string = flag.String("parentOne", "", "Name of first parental genome``")
	var parentTwo *string = flag.String("parentTwo", "", "Name of second parental genome``")
	//TODO: Look into using the groups struct in the popgen package
	var list *string = flag.String("byname", "", "Filter samples of interest by providing a name`.txt` file containing a list of sample names, one name per line")

	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		if *sampleName && len(flag.Args()) == 1 {
			file := fileio.EasyOpen(flag.Arg(0))
			defer file.Close()
			header := vcf.ReadHeader(file)
			fmt.Printf("%s", vcf.PrintSampleNames(header))
		} else {
			flag.Usage()
			fmt.Printf("\nExamples:\nAllele Specific Filter:\n./filterGenotypes -f1 name -parentOne name -parentTwo name input.vcf output.vcf\n\nView sample names:\n./filterGenotypes -samples file.vcf\n\n")
			log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
		}
	} else {
		input, output := flag.Arg(0), flag.Arg(1)
		if strings.HasSuffix(*list, ".txt") {
			samples := fileio.Read(*list)

			gvcf := vcf.GoReadGVcf(input)
			writer := fileio.EasyCreate(output)

			//write header
			vcf.WriteMultiSamplesHeader(writer, gvcf.Header, samples)
			vcf.ByNames(gvcf, samples, writer)

			writer.Close()

		} else if *parentOne == "" || *parentTwo == "" || *f1Genome == "" {
			log.Fatalf("Error: Must provide exactly 2 parents and 1 F1 sample...\n")
		} else {
			reader, header := vcf.GoReadToChan(input)
			sampleHash := vcf.HeaderToMaps(header)

			var parentalOne, parentalTwo, fOne int16 = sampleHash.GIndex[*parentOne], sampleHash.GIndex[*parentTwo], sampleHash.GIndex[*f1Genome]
			writer := fileio.EasyCreate(output)
			vcf.NewWriteHeader(writer, header)


			defer writer.Close()
			log.SetFlags(0)
			for each := range reader {
				if vcf.ASFilter(each, parentalOne, parentalTwo, fOne) {
					//TODO: Need to chanch logic for this log print feature.
					//vcf.PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
					vcf.WriteVcf(writer, each)
				}
			}
		}
	}
}
