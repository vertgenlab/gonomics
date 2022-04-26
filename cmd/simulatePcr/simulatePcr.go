// Command Group: "Data Simulation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
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

func simulatePcr(fwdPrimers primerSeqs, ref string, out string, maxLen int) {
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
		go findAmplicons(chrom, seq, fwdPrimers, revPrimers, maxLen, ans, wg) // spawn goroutine for each chrom
	}

	go func() {
		wg.Wait()
		close(ans)
	}()

	var recordsWritten int
	o := fileio.EasyCreate(out)
	for b := range ans {
		bed.WriteBed(o, b)
		recordsWritten++
	}

	err := o.Close()
	exception.PanicOnErr(err)

	log.Printf("found %d potential products", recordsWritten) // prefer an explicit report rather than empty file for zero amplicons found
}

// findAmplicons determins the amplicons resulting from the input primers for a single given template sequence
func findAmplicons(chrom, template string, fwdPrimers, revPrimers []string, maxLen int, ans chan<- bed.Bed, chromWg *sync.WaitGroup) {
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

	calcProducts(chrom, fwdSites, revSites, maxLen, ans)
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
		if fwd { // the start for forward primers is the base after the 3' end of the primer
			currPos += len(primer)
		}
		send <- primingSite{primer: primer, start: offset + currPos}
		if !fwd { // add to the rev pos after sending answer to be consistent with fwd
			currPos += len(primer)
		}

		offset += currPos
		template = template[currPos:]
	}
	wg.Done()
}

// calcProducts inputs forward and reverse priming sites and sends the resulting products to the ans channel.
func calcProducts(chrom string, fwdSites, revSites []primingSite, maxLen int, ans chan<- bed.Bed) {
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
	out := flag.String("o", "stdout", "Output bed file location.")
	maxLen := flag.Int("maxProductSize", 1000, "Limits the size of output products. Set to zero for no limit.")
	maxProcs := flag.Int("maxCPU", 0, "Limits the number of CPUs available to the program. Set to 0 for no limit.")
	flag.Parse()
	flag.Usage = usage

	runtime.GOMAXPROCS(*maxProcs)

	if len(primers) == 0 {
		usage()
		log.Fatal("ERROR: Primer sequence must be declared at least once using -p")
	}

	if *ref == "" {
		usage()
		log.Fatal("ERROR: Must declare template sequence with -t")
	}

	if *maxLen == 0 {
		*maxLen = math.MaxInt
	}

	simulatePcr(primers, *ref, *out, *maxLen)
}
