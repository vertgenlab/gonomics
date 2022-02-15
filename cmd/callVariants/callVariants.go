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
	"sync"
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

func callVariants(experimentalFiles, normalFiles []string, refFile, outFile string, maxP, minAf, maxAf, maxStrandBias float64, minCoverage, minMapQ, minAltReads, workers int) {
	out := fileio.EasyCreate(outFile)
	outHeader := makeOutputHeader(append(experimentalFiles, normalFiles...))
	vcf.NewWriteHeader(out, outHeader)

	readFilters := getReadFilters(uint8(minMapQ))
	pileFilters := getPileFilters(minCoverage)
	expHeaders, expPiles := startPileup(experimentalFiles, readFilters, pileFilters)
	normHeaders, normPiles := startPileup(normalFiles, readFilters, pileFilters)

	isExperimental := make([]bool, len(experimentalFiles)+len(normalFiles))
	for i := 0; i < len(experimentalFiles); i++ {
		isExperimental[i] = true
	}

	err := checkHeadersMatch(append(expHeaders, normHeaders...))
	exception.FatalOnErr(err)

	synced := sam.GoSyncPileups(append(expPiles, normPiles...)...)

	// start goroutines to process synced piles independently and send results to writeChan
	writeChan := make(chan vcf.Vcf, 1000)
	wg := new(sync.WaitGroup)
	for i := 0; i < workers; i++ {
		wg.Add(1)
		go startWorker(wg, writeChan, synced, len(experimentalFiles), expHeaders[0], fasta.NewSeeker(refFile, ""), maxP, minAf, maxAf, maxStrandBias, minCoverage, minAltReads)
	}

	// spawn a goroutine that waits until all the goroutines are done then closes writeChan.
	go func(*sync.WaitGroup, chan vcf.Vcf) {
		wg.Wait()
		close(writeChan)
	}(wg, writeChan)

	// write the vcf files as they come in
	for v := range writeChan {
		vcf.WriteVcf(out, v)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func startWorker(wg *sync.WaitGroup, writeChan chan<- vcf.Vcf, synced <-chan []sam.Pile, numExpFiles int, header sam.Header, ref *fasta.Seeker, maxP, minAf, maxAf, maxStrandBias float64, minCoverage, minAltReads int) {
	var keepVar bool
	for piles := range synced {
		var v vcf.Vcf
		v, keepVar = getVariant(piles[:numExpFiles], piles[numExpFiles:], header, ref, maxP, minAf, maxAf, maxStrandBias, minCoverage, minAltReads)
		if keepVar {
			writeChan <- v
		}
	}
	err := ref.Close()
	exception.PanicOnErr(err)
	wg.Done()
}

// startPileup for each input file
func startPileup(files []string, readFilters []func(s sam.Sam) bool, pileFilters []func(p sam.Pile) bool) (headers []sam.Header, piles []<-chan sam.Pile) {
	headers = make([]sam.Header, len(files))
	piles = make([]<-chan sam.Pile, len(files))

	var samChan <-chan sam.Sam
	for i := range files {
		samChan, headers[i] = sam.GoReadToChan(files[i])
		if len(headers[i].Text) == 0 {
			log.Fatal("ERROR: sam/bam files must have headers")
		}
		piles[i] = sam.GoPileup(samChan, headers[i], false, readFilters, pileFilters)
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

func getPileFilters(minCoverage int) []func(p sam.Pile) bool {
	var answer []func(p sam.Pile) bool
	answer = append(answer, func(p sam.Pile) bool {
		return calcDepth(p) >= minCoverage
	})
	return answer
}

func getReadFilters(minMapQ uint8) []func(s sam.Sam) bool {
	var answer []func(s sam.Sam) bool
	answer = append(answer, func(s sam.Sam) bool {
		return s.MapQ >= minMapQ
	})
	return answer
}

func main() {
	var experimentalFiles, normalFiles inputFiles
	flag.Var(&experimentalFiles, "i", "Input experimental files. May be declared more than once (.bam, .sam)")
	flag.Var(&normalFiles, "n", "Input normal files. May be declared more than once (.bam, .sam)")

	maxP := flag.Float64("p", 0.001, "Maximum p-value for output")
	minAf := flag.Float64("minAF", 0.01, "Minimum allele frequency of variants")
	maxAf := flag.Float64("maxAF", 1, "Maximum allele frequency of variants")
	maxStrandBias := flag.Float64("maxStrandBias", 0.9, "Maximum allowable strand imbalance of reads supporting the alternate allele. "+
		"When set to 0.9, any variants where >90% of reads supporting the alternate allele are on the forward/reverse strand, that variant will be discarded.")
	minCoverage := flag.Int("minCoverage", 10, "Minimum coverage for site to be considered.")
	minMapQ := flag.Int("minMapQ", 10, "Minimum mapping quality for read to be considered")
	minAltReads := flag.Int("minAltReads", 1, "Minimum reads supporting the alternate allele for variant to be considered")
	reffile := flag.String("r", "", "Reference fasta file (.fa). Must be indexed (.fai)")
	outfile := flag.String("o", "stdout", "Output file (.vcf)")
	workers := flag.Int("t", 1, "Number of threads used for calling variants. Note that additional threads are automatically generated for io.")
	flag.Parse()
	flag.Usage = usage

	if len(experimentalFiles) == 0 {
		flag.Usage()
		log.Fatalln("ERROR: must declare at least 1 experimental sample with -i")
	}

	if *reffile == "" {
		flag.Usage()
		log.Fatalln("ERROR: must input an indexed fasta file with -r")
	}

	callVariants(experimentalFiles, normalFiles, *reffile, *outfile, *maxP, *minAf, *maxAf, *maxStrandBias, *minCoverage, *minMapQ, *minAltReads, *workers)
}
