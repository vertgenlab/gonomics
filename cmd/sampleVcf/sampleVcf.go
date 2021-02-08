package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math/rand"
)

//sampleVcf takes a VCF file and returns a random subset of variants to an output VCF file. Can also retain a random subset of alleles from gVCF data (diploid, does not break allele pairs)
func sampleVcf(inFile string, outFile string, numVariants int, numSamples int, randSeed bool, setSeed int64) {
	common.RngSeed(randSeed, setSeed)
	records, header := vcf.ReadWithHeader(inFile)

	var sampleList []string = vcf.HeaderGetSampleList(header)

	if numVariants > len(records) {
		log.Fatalf("The number of requested sampled variants is greater than the number of variants in the input file.")
	}

	//Shuffle the vcf records, our subset will be composed to the first entries in the shuffled order.
	rand.Shuffle(len(records), func(i, j int) { records[i], records[j] = records[j], records[i] })
	//DEBUG:fmt.Printf("lenRecords before slice: %v.\n", len(records))
	records = records[:numVariants] //keep only as many results as specified.
	//DEBUG: fmt.Printf("lenRecords after slice: %v.\n", len(records))

	if numSamples > 0 {
		if numSamples > len(records[0].Samples) {
			log.Fatalf("More samples were requested than were present in the input VCF file.")
		}
		var sequentialSlice []int = getSequentialSlice(len(records[0].Samples))
		rand.Shuffle(len(sequentialSlice), func(i, j int) { sequentialSlice[i], sequentialSlice[j] = sequentialSlice[j], sequentialSlice[i] })
		sequentialSlice = sequentialSlice[:numSamples] //now we have a list of samples to keep from each variant.

		var outHeaderSampleList []string = make([]string, 0)
		for _, i := range sequentialSlice {
			outHeaderSampleList = append(outHeaderSampleList, sampleList[i])
		}

		vcf.HeaderUpdateSampleList(header, outHeaderSampleList)

		var outSamples []vcf.GenomeSample

		for _, i := range records {
			outSamples = make([]vcf.GenomeSample, 0)
			for _, j := range sequentialSlice {
				outSamples = append(outSamples, i.Samples[j])
			}
			i.Samples = outSamples
		}
	}

	out := fileio.EasyCreate(outFile)
	defer out.Close()
	vcf.NewWriteHeader(out, header)
	vcf.WriteVcfToFileHandle(out, records)
}

//returns a slice where the value is the index. Answer is of length n. ex (4) returns [0 1 2 3]
func getSequentialSlice(n int) []int {
	var answer []int = make([]int, n)
	for i := 0; i < n; i++ {
		answer[i] = i
	}
	return answer
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
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
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

	sampleVcf(inFile, outFile, *numVariants, *numSamples, *randSeed, *setSeed)
}
