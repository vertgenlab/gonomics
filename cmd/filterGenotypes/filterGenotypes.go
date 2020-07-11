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
		"filterGenotypes - a gonomics tool to filter genotyped Vcfs containing multiple samples\n\n" +
			"Usage:\n" +
			"  ./filterGenotypes [options] input.vcf output.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	var parental arrayFlags
	flag.Var(&parental, "parent", "Define of two parential that appears homozygous in the genotype Vcf. Must match Vcf Header exactly and include second `name -parent two -f1 name`")
	var f1Genome *string = flag.String("f1", "", "F1 hybrid sample that appears heterozygous in genotype Vcf. Must match Vcf Header exactly and include `name -parent one -parent two`")
	var sampleName *bool = flag.Bool("samples", false, "Get names of samples that appear in Vcf header. (Default: /dev/stdout)")

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
			log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
		}
	} else {
		input, output := flag.Arg(0), flag.Arg(1)
		if strings.HasSuffix(*list, ".txt") {
			samples := fileio.Read(*list)

			gvcf := vcf.ReadGVcf(input)
			writer := fileio.EasyCreate(output)

			//write header
			vcf.WriteMultiSamplesHeader(writer, gvcf.Header, samples)
			vcf.ByNames(gvcf, samples, writer)

			writer.Close()

		} else if len(parental) != 2 {
			log.Fatalf("Error: Must provide exactly 2 parents, found %d...\n", len(parental))
		} else {
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
