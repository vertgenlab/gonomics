// Command Group: "Variant Calling & Annotation"

package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"path/filepath"
	"strings"
	"time"
)

func usage() {
	fmt.Print(
		"callVariants - A tool to find variation between multiple alignment files.\n\n" +
			"Usage:\n" +
			"  callVariants [options] -i file1.bam -i file2.bam -n normal.bam -r reference.fasta\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func callVariants(experimentalFiles, normalFiles []string, reffile, outfile string, maxP, minAf float64, minCoverage int) {
	out := fileio.EasyCreate(outfile)
	outHeader := makeOutputHeader(append(experimentalFiles, normalFiles...))
	vcf.NewWriteHeader(out, outHeader)
	ref := fasta.NewSeeker(reffile, "")

	expHeaders, expPiles := startPileup(experimentalFiles)
	normHeaders, normPiles := startPileup(normalFiles)

	isExperimental := make([]bool, len(experimentalFiles)+len(normalFiles))
	for i := 0; i < len(experimentalFiles); i++ {
		isExperimental[i] = true
	}

	err := checkHeadersMatch(append(expHeaders, normHeaders...))
	exception.FatalOnErr(err)

	synced := sam.GoSyncPileups(append(expPiles, normPiles...)...)

	var v vcf.Vcf
	var keepVar bool
	for piles := range synced {
		v, keepVar = getVariant(piles[:len(experimentalFiles)], piles[len(experimentalFiles):], ref, expHeaders[0].Chroms, maxP, minAf, minCoverage)
		if keepVar {
			vcf.WriteVcf(out, v)
		}
	}

	err = ref.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

// startPileup for each input file
func startPileup(files []string) (headers []sam.Header, piles []<-chan sam.Pile) {
	headers = make([]sam.Header, len(files))
	piles = make([]<-chan sam.Pile, len(files))

	var samChan <-chan sam.Sam
	for i := range files {
		samChan, headers[i] = sam.GoReadToChan(files[i])
		if len(headers[i].Text) == 0 {
			log.Fatal("ERROR: sam/bam files must have headers")
		}
		piles[i] = sam.GoPileup(samChan, headers[i], false, nil, nil)
	}
	return
}

// checkHeadersMatch verifies that all input files use the same reference
func checkHeadersMatch(headers []sam.Header) error {
	ref := headers[0].Chroms
	for i := 1; i < len(headers); i++ {
		if len(headers[i].Chroms) != len(ref) {
			return errors.New("ERROR: reference chromosomes in input files must match")
		}
		for j := range headers[i].Chroms {
			if headers[i].Chroms[j] != ref[j] {
				return errors.New("ERROR: reference chromosomes in input files must match and be in the same order")
			}
		}
	}
	return nil
}

// makeOutputHeader produces a header for the output vcf file
func makeOutputHeader(filenames []string) vcf.Header {
	var header vcf.Header
	sampleNames := make([]string, len(filenames))
	for i := range filenames {
		sampleNames[i] = strings.TrimSuffix(filepath.Base(filenames[i]), filepath.Ext(filenames[i]))
	}
	t := time.Now()
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, "##fileDate="+t.Format("20060102"))
	header.Text = append(header.Text, "##source=github.com/vertgenlab/gonomics")
	header.Text = append(header.Text, "##phasing=none")
	header.Text = append(header.Text, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	header.Text = append(header.Text, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">")
	header.Text = append(header.Text, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth of Each Allele\">")
	header.Text = append(header.Text, "##FORMAT=<ID=PV,Number=A,Type=Floatg,Description=\"p value for Each Alternate Allele\">")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", strings.Join(sampleNames, "\t")))
	return header
}

// inputFiles is a custom type that gets filled by flag.Parse()
type inputFiles []string

// String to satisfy flag.Value interface
func (i *inputFiles) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *inputFiles) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var experimentalFiles, normalFiles inputFiles
	flag.Var(&experimentalFiles, "i", "Input experimental files. May be declared more than once (.bam, .sam)")
	flag.Var(&normalFiles, "n", "Input normal files. May be declared more than once (.bam, .sam)")

	maxP := flag.Float64("p", 0.001, "Maximum p-value for output")
	minAf := flag.Float64("af", 0.01, "Minimum allele frequency of variants")
	minCoverage := flag.Int("minCoverage", 10, "Minimum coverage for site to be considered.")
	reffile := flag.String("r", "", "Reference fasta file (.fa). Must be indexed (.fai)")
	outfile := flag.String("o", "stdout", "Output file (.vcf)")
	flag.Parse()
	flag.Usage = usage

	if len(experimentalFiles) == 0 {
		flag.Usage()
		log.Fatalln("ERROR: must declare at least 1 experimental sample with -i")
	}

	callVariants(experimentalFiles, normalFiles, *reffile, *outfile, *maxP, *minAf, *minCoverage)
}
