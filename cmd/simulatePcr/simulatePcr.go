// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"math"
	"runtime"
	"sort"
	"strings"
	"sync"
)

func usage() {
	fmt.Print(
		"simulatePcr - predict amplicon products generated from PCR.\n\n" +
			"Usage:\n" +
			"  simulatePcr [options] -r reference.fasta -p FOR_PRIMER_SEQ -p REV_PRIMER_SEQ\n" +
			"Note that simulatePcr assumes only perfect matches will amplify and only the smallest of" +
			"competing overlapping amplicons will be reported\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

// primerSeqs is a custom type that gets filled by flag.Parse()
type primerSeqs []string

// String to satisfy flag.Value interface
func (i *primerSeqs) String() string {
	return strings.Join(*i, " ")
}

// Set to satisfy flag.Value interface
func (i *primerSeqs) Set(value string) error {
	*i = append(*i, value)
	return nil
}

type primingSite struct {
	primer string
	//fwd    bool // true if the primer seq matches the + strand (i.e. a forward primer)
	start int // first base of the primer for reverse strand, first base AFTER the primer for forward strand
	// such that the product is seq[fwd.start:rev.start]
}

func simulatePcr(fwdPrimers primerSeqs, ref string, outBed, outFastq string, maxLen int, includePrimer bool) {
	// we want strings here to use the strings.Index function
	var template map[string]string = fasta.ReadToString(ref)

	// add the reverse complement of each primer to the slice of primers
	revPrimers := make([]string, len(fwdPrimers))
	var complement []dna.Base
	for i := 0; i < len(fwdPrimers); i++ {
		complement = dna.StringToBases(fwdPrimers[i])
		dna.ReverseComplement(complement)
		revPrimers[i] = dna.BasesToString(complement)
	}

	ans := make(chan bed.Bed, 1000)
	wg := new(sync.WaitGroup)
	for chrom, seq := range template {
		wg.Add(1)
		go findAmplicons(chrom, seq, fwdPrimers, revPrimers, maxLen, includePrimer, ans, wg) // spawn goroutine for each chrom
	}

	go func() {
		wg.Wait()
		close(ans)
	}()

	fqQualSlice := make([]uint8, maxLen)
	for i := range fqQualSlice {
		fqQualSlice[i] = 40 // set to perfect score
	}
	var recordsWritten int
	var outBedHandle, outFqHandle *fileio.EasyWriter
	if outBed != "" {
		outBedHandle = fileio.EasyCreate(outBed)
	}
	if outFastq != "" {
		outFqHandle = fileio.EasyCreate(outFastq)
	}
	var currFq fastq.Fastq
	for b := range ans {
		if outBed != "" {
			bed.WriteBed(outBedHandle, b)
		}
		if outFastq != "" {
			currFq.Name = fmt.Sprintf("%s:%d-%d_%s", b.Chrom, b.ChromStart, b.ChromEnd, b.Name)
			currFq.Seq = dna.StringToBases(template[b.Chrom][b.ChromStart:b.ChromEnd])
			currFq.Qual = fqQualSlice[:len(currFq.Seq)]
			fastq.WriteToFileHandle(outFqHandle, currFq)
		}
		recordsWritten++
	}

	var err error
	if outBed != "" {
		err = outBedHandle.Close()
		exception.PanicOnErr(err)
	}
	if outFastq != "" {
		err = outFqHandle.Close()
		exception.PanicOnErr(err)
	}
	log.Printf("found %d potential products", recordsWritten) // prefer an explicit report rather than empty file for zero amplicons found
}

// findAmplicons determines the amplicons resulting from the input primers for a single given template sequence
func findAmplicons(chrom, template string, fwdPrimers, revPrimers []string, maxLen int, includePrimer bool, ans chan<- bed.Bed, chromWg *sync.WaitGroup) {
	template = strings.ToUpper(template)
	var fwdSites, revSites []primingSite
	fwdChan := make(chan primingSite, 100)
	revChan := make(chan primingSite, 100)

	primerWg := new(sync.WaitGroup)
	for i := range fwdPrimers { // spawn a goroutine for each primer
		primerWg.Add(2)
		go findPrimingSites(template, fwdPrimers[i], true, fwdChan, primerWg)
		go findPrimingSites(template, revPrimers[i], false, revChan, primerWg)
	}

	compWg := new(sync.WaitGroup)
	compWg.Add(2)
	go func() { // compile fwd sites
		for s := range fwdChan {
			fwdSites = append(fwdSites, s)
		}
		compWg.Done()
	}()

	go func() { // compile rev sites
		for s := range revChan {
			revSites = append(revSites, s)
		}
		compWg.Done()
	}()

	primerWg.Wait()
	close(fwdChan)
	close(revChan)
	compWg.Wait()

	sort.Slice(fwdSites, func(i, j int) bool {
		return fwdSites[i].start < fwdSites[j].start
	})

	sort.Slice(revSites, func(i, j int) bool {
		return revSites[i].start < revSites[j].start
	})

	calcProducts(chrom, fwdSites, revSites, maxLen, includePrimer, ans)
	chromWg.Done()
}

// findPrimingSites determines where the input primer will bind to the input template.
func findPrimingSites(template, primer string, fwd bool, send chan<- primingSite, wg *sync.WaitGroup) {
	var currPos, offset int
	for len(template) > 0 {
		currPos = strings.Index(template, primer)
		if currPos == -1 { // no priming site found
			break
		}

		if fwd { // fwd first base is last base after the primer
			currPos += len(primer)
		}

		send <- primingSite{primer: primer, start: offset + currPos}

		if !fwd { // increment the reverse primer after sending to match fwd primer
			currPos += len(primer)
		}

		offset += currPos
		template = template[currPos:]
	}
	wg.Done()
}

// calcProducts inputs forward and reverse priming sites and sends the resulting products to the ans channel.
func calcProducts(chrom string, fwdSites, revSites []primingSite, maxLen int, includePrimer bool, ans chan<- bed.Bed) {
	if len(fwdSites) == 0 || len(revSites) == 0 {
		return
	}

	var revIdx, start, end int
	var revBases []dna.Base
	for i := 0; i < len(fwdSites); i++ {
		for fwdSites[i].start >= revSites[revIdx].start {
			revIdx++ // fwd site is after rev site, increment rev until after fwd
			if revIdx >= len(revSites) {
				return
			}
		}

		if i+1 < len(fwdSites) && fwdSites[i+1].start < revSites[revIdx].start {
			continue // next fwd site is closer in, move to next site
		}

		start = fwdSites[i].start
		end = revSites[revIdx].start

		if includePrimer {
			start -= len(fwdSites[i].primer)
			end += len(revSites[revIdx].primer)
		}
		if end-start <= maxLen {
			revBases = dna.StringToBases(revSites[revIdx].primer)
			dna.ReverseComplement(revBases) // give input primer sequence in report

			ans <- bed.Bed{
				Chrom:             chrom,
				ChromStart:        start,
				ChromEnd:          end,
				Name:              fwdSites[i].primer + "+" + dna.BasesToString(revBases),
				FieldsInitialized: 4,
			}
		}
	}
}

func main() {
	var primers primerSeqs
	flag.Var(&primers, "p", "Primer sequence (5'->3', ACGT only). Must be declared at least once.")
	ref := flag.String("t", "", "Template sequence to use for amplification (fasta).")
	outBed := flag.String("bed", "", "Output bed file location.")
	outFastq := flag.String("fastq", "", "Output fastq file file location.")
	maxLen := flag.Int("maxProductSize", 1000, "Limits the size of output products. Set to zero for no limit.")
	maxProcs := flag.Int("maxCPU", 0, "Limits the number of CPUs available to the program. Set to 0 for no limit.")
	includePrimer := flag.Bool("includePrimer", true, "Include primer sequence in output bed/fastq.")
	flag.Parse()
	flag.Usage = usage

	runtime.GOMAXPROCS(*maxProcs)

	if len(primers) == 0 {
		usage()
		log.Fatal("ERROR: primer sequence must be declared at least once using -p")
	}

	if *ref == "" {
		usage()
		log.Fatal("ERROR: must declare template sequence with -t")
	}

	if *maxLen == 0 {
		*maxLen = math.MaxInt
	}

	if *outBed == "" && *outFastq == "" {
		usage()
		log.Fatal("ERROR: must declare at least one output file (-bed or -fastq)")
	}

	simulatePcr(primers, *ref, *outBed, *outFastq, *maxLen, *includePrimer)
}
