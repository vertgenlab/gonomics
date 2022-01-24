// Command Group: "Variant Calling & Annotation"

package main

import (
	"errors"
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

func callVariants(experimentalFiles, normalFiles []string, reffile, outfile string, maxP float64) {
	out := fileio.EasyCreate(outfile)
	outHeader := vcf.NewHeader(strings.Join(append(experimentalFiles, normalFiles...), "\t"))
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

	synced := sam.GoSyncPileups(append(expPiles, normPiles...))

	var v vcf.Vcf
	var keepVar bool
	for piles := range synced {
		v, keepVar = getVariant(piles[:len(experimentalFiles)], piles[len(experimentalFiles):], ref, expHeaders[0].Chroms, maxP)
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
	flag.Var(&experimentalFiles, "i", "Input experimental files , may be declared more than once (.bam, .sam)")
	flag.Var(&normalFiles, "n", "Input normal files, may be declared more than once (.bam, .sam)")

	maxP := flag.Float64("p", 0.05, "Maximum p-value for output")
	reffile := flag.String("r", "", "Reference fasta file (.fa). Must be indexed (.fai)")
	outfile := flag.String("o", "stdout", "Output file (.vcf)")
	flag.Parse()

	callVariants(experimentalFiles, normalFiles, *reffile, *outfile, *maxP)
}
