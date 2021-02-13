package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math"
	"strings"
)

func vcfFilter(infile string, outfile string, groupFile string, chrom string, minPos int, maxPos int, ref string, alt string, minQual float64) {
	ch, header := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	defer out.Close()
	altSlice := strings.Split(alt, ",")

	var samplesToKeep []int = make([]int, 0) //this var holds all of the indices from samples (defined below as the sample list in the header) that we want to keep in the output file.

	if groupFile != "" {
		groups := popgen.ReadGroups(groupFile)
		samples := vcf.HeaderGetSampleList(header)

		for i := 0; i < len(samples); i++ {
			if popgen.GroupsContains(groups, samples[i]) {
				samplesToKeep = append(samplesToKeep, i)
			}
		}
		outSamples := filterHeaderSamplesToKeep(samples, samplesToKeep)
		vcf.HeaderUpdateSampleList(header, outSamples)
	}

	vcf.NewWriteHeader(out.File, header)

	for v := range ch {
		if vcf.Filter(v, chrom, minPos, maxPos, ref, altSlice, minQual) {
			v.Samples = filterRecordsSamplesToKeep(v.Samples, samplesToKeep)
			vcf.WriteVcf(out.File, v)
		}
	}
}

func filterRecordsSamplesToKeep(recordSamples []vcf.GenomeSample, samplesToKeep []int) []vcf.GenomeSample {
	var answer []vcf.GenomeSample
	for i := 0; i < len(samplesToKeep); i++ {
		answer = append(answer, recordSamples[samplesToKeep[i]])
	}

	return answer
}

func filterHeaderSamplesToKeep(samples []string, samplesToKeep []int) []string {
	var answer []string
	for i := 0; i < len(samplesToKeep); i++ {
		answer = append(answer, samples[samplesToKeep[i]])
	}
	return answer
}

func usage() {
	fmt.Print(
		"vcfFilter\n" +
			"Usage:\n" +
			"vcfFilter input.vcf output.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var groupFile *string = flag.String("groupFile", "", "Retains alleles from individuals in the input group file.")
	var minPos *int = flag.Int("minPos", math.MinInt64, "Specifies the minimum position of the variant.")
	var maxPos *int = flag.Int("maxPos", math.MaxInt64, "Specifies the maximum position of the variant.")
	var minQual *float64 = flag.Float64("minQual", 0.0, "Specifies the minimum quality score.")
	var ref *string = flag.String("ref", "", "Specifies the reference field.")
	var alt *string = flag.String("alt", "", "Specifies the alt field.")

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

	vcfFilter(infile, outfile, *groupFile, *chrom, *minPos, *maxPos, *ref, *alt, *minQual)
}
