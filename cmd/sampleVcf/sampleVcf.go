// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

//sampleVcf takes a VCF file and returns a random subset of variants to an output VCF file. Can also retain a random subset of alleles from gVCF data (diploid, does not break allele pairs)
func sampleVcf(inFile string, outFile string, numVariants int, numSamples int, setSeed int64) {
	rand.Seed(setSeed)
	records, header := vcf.Read(inFile)
	records, header = vcf.SampleVcf(records, header, numVariants, numSamples)

	out := fileio.EasyCreate(outFile)
	vcf.NewWriteHeader(out, header)
	vcf.WriteVcfToFileHandle(out, records)
	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"sampleVcf - Returns a sample from a VCF file with a specified number of results.\n" +
			"Usage:\n" +
			"sampleVcf input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var numVariants *int = flag.Int("numVariants", 1, "Specifies the number of variants to return from file.")
	var numSamples *int = flag.Int("numSamples", -1, "Specifies the number of samples (diploid allele pairs) to retain in the output sample.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	sampleVcf(inFile, outFile, *numVariants, *numSamples, *setSeed)
}
