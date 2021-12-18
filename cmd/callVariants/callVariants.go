// Command Group: "Variant Calling & Annotation"

package main

import (
	"errors"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"strings"
)

func callVariants(experimentalFiles, normalFiles []string, outfile string) {
	expHeaders, expPiles := startPileup(experimentalFiles)
	normHeaders, normPiles := startPileup(normalFiles)
	isExperimental := make([]bool, len(experimentalFiles)+len(normalFiles))
	for i := 0; i < len(experimentalFiles); i++ {
		isExperimental[i] = true
	}

	err := checkHeadersMatch(append(expHeaders, normHeaders...))
	exception.FatalOnErr(err)

	synced := sam.GoSyncPileups(append(expPiles, normPiles...))
	for piles := range synced {
		fmt.Println(piles)
	}
}

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

type inputFiles []string

func (i *inputFiles) String() string {
	return strings.Join(*i, " ")
}

func (i *inputFiles) Set(value string) error {
	*i = append(*i, value)
	return nil
}

func main() {
	var experimentalFiles, normalFiles inputFiles
	flag.Var(&experimentalFiles, "i", "Input experimental files , may be declared more than once (.bam, .sam)")
	flag.Var(&normalFiles, "n", "Input normal files, may be declared more than once (.bam, .sam)")

	outfile := flag.String("o", "stdout", "Output file (.vcf)")
	flag.Parse()

	callVariants(experimentalFiles, normalFiles, *outfile)
}
