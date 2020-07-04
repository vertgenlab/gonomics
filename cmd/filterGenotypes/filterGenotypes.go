package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func usage() {
	fmt.Print(
		"goVcf - a gonomics tool to filter genotyped Vcfs containing multiple samples\n\n" +
			"Usage:\n" +
			"  ./goVcf [options] input output [.vcf/.gz] \n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	var parental arrayFlags
	flag.Var(&parental, "parent", "Define of two parential that appears homozygous in the genotype Vcf\nMust match Vcf Header exactly``")
	var f1Genome *string = flag.String("f1", "", "Define of f1 hybrid sample that appears heterozygous in the genotype Vcf\nMust match Vcf Header exactly``")
	var sampleName *bool = flag.Bool("sampleNames", false, "Get names of samples that appear in Vcf header\nDoes not require output file")

	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		if *sampleName && len(flag.Args()) == 1 {
			file := fileio.EasyOpen(flag.Arg(0))
			defer file.Close()
			header := vcf.ReadHeader(file)
			log.Printf("\n%s", vcf.PrintSampleNames(header))
		} else {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
		}
	} else {
		input, output := flag.Arg(0), flag.Arg(1)
		if len(parental) != 2 {
			log.Fatalf("Error: Must provide exactly 2 parents, found %d...\n", len(parental))
		}
		file := fileio.EasyOpen(input)
		defer file.Close()

		header := vcf.ReadHeader(file)
		sampleHash := vcf.HeaderToMaps(header)

		var parentalOne, parentalTwo, fOne int16 = sampleHash.IndexAllele[parental[0]], sampleHash.IndexAllele[parental[1]], sampleHash.IndexAllele[*f1Genome]
		writer := fileio.EasyCreate(output)
		vcf.NewWriteHeader(writer, header)
		reader := make(chan *vcf.Vcf)
		go vcf.ReadToChan(file, reader)

		defer writer.Close()
		log.SetFlags(0)
		for each := range reader {
			if vcf.ASFilter(each, parentalOne, parentalTwo, fOne) {
				vcf.PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
				vcf.WriteVcf(writer, each)
			}
		}
	}
}

//Define flag value as an array. Used to define parent genomes.
type arrayFlags []string

func (i *arrayFlags) String() string {
	return "my string representation"
}
func (i *arrayFlags) Set(value string) error {
	*i = append(*i, value)
	return nil
}
