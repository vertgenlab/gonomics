// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/popgen"
	"github.com/vertgenlab/gonomics/vcf"
)

// coordinate is a chrom/pos pair, referring to a specific place in the genome.
type coordinate struct {
	Chr string
	Pos int
}

func getSitesSeen(filename string) map[coordinate]uint8 {
	var sitesSeen map[coordinate]uint8 = make(map[coordinate]uint8, 0)
	var records <-chan vcf.Vcf
	records, _ = vcf.GoReadToChan(filename)
	var currentCoord coordinate = coordinate{"", 0}
	for v := range records {
		currentCoord.Chr = v.Chr
		currentCoord.Pos = v.Pos
		sitesSeen[currentCoord]++
	}
	return sitesSeen
}

func rmClusteredRecords(input <-chan vcf.Vcf, output chan<- vcf.Vcf, minDist int, totalChan, removedChan chan<- int) {
	var prev vcf.Vcf
	var canSend bool
	var total, removed int
	for v := range input {
		total++
		if prev.Chr == "" {
			prev = v
			canSend = true
			continue
		}

		if v.Pos < prev.Pos && v.Chr == prev.Chr {
			log.Fatalf("ERROR: input vcf is not sorted. Offending records:\n%s\n%s", prev, v)
		}

		// if onto a new chrom, send previous if applicable, set v as prev, and move on to next record
		if v.Chr != prev.Chr {
			if canSend {
				output <- prev
			} else {
				removed++
			}
			canSend = true
			prev = v
			continue
		}

		if v.Pos-prev.Pos < minDist { // too close, don't send
			canSend = false
			prev = v
			removed++
			continue
		}

		// if we reach this point prev is clear to send and v is ok to send in the next loop if the subsequent variant is not too close
		if canSend {
			output <- prev
		} else {
			removed++
		}
		prev = v
		canSend = true
	}

	// send final record if passing
	if canSend {
		output <- prev
	} else {
		removed++
	}

	totalChan <- total
	removedChan <- removed

	close(totalChan)
	close(removedChan)
	close(output)
}

func vcfFilter(infile string, outfile string, c criteria, groupFile string, parseFormat bool, parseInfo bool, setSeed int64) (total, removed int) {
	rand.Seed(setSeed)
	var records <-chan vcf.Vcf
	var totalChan, removedChan chan int
	var header vcf.Header
	var err error
	var currentCoord coordinate
	var sitesSeen map[coordinate]uint8 = make(map[coordinate]uint8, 0) //uint8 is the number of times this site is seen in the vcf file.

	if c.biAllelicOnly {
		sitesSeen = getSitesSeen(infile)
	}

	records, header = vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	tests := getTests(c, header)

	if c.minDist > 0 {
		// hijack input channel so only sites passing the minDist filter are passed to the rest of the function
		passedRecords := make(chan vcf.Vcf, 100)
		totalChan = make(chan int, 1)   // to get accurate total count
		removedChan = make(chan int, 1) // to get accurate removed count
		go rmClusteredRecords(records, passedRecords, c.minDist, totalChan, removedChan)
		records = passedRecords
	}

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
	vcf.NewWriteHeader(out, header)

	for v := range records {
		total++
		if groupFile != "" {
			v.Samples = filterRecordsSamplesToKeep(v.Samples, samplesToKeep)
		}

		if parseFormat {
			v = vcf.ParseFormat(v, header)
		}

		if c.biAllelicOnly {
			currentCoord = coordinate{v.Chr, v.Pos}
			if sitesSeen[currentCoord] < 1 {
				log.Panicf("Current variant not found in sitesSeen map. Something went horribly wrong. %v.\n", v)
			} else if sitesSeen[currentCoord] > 1 {
				removed++
				continue
			}
		}

		if parseInfo {
			v = vcf.ParseInfo(v, header)
		}

		if !passesTests(v, tests) {
			removed++
			continue
		}

		vcf.WriteVcf(out, v)
	}

	if c.minDist > 0 {
		total = <-totalChan
		removed += <-removedChan
	}

	err = out.Close()
	exception.PanicOnErr(err)
	return
}

func filterRecordsSamplesToKeep(recordSamples []vcf.Sample, samplesToKeep []int) []vcf.Sample {
	var answer []vcf.Sample
	for _, v := range samplesToKeep {
		answer = append(answer, recordSamples[v])
	}

	return answer
}

func filterHeaderSamplesToKeep(samples []string, samplesToKeep []int) []string {
	var answer []string
	for _, v := range samplesToKeep {
		answer = append(answer, samples[v])
	}

	return answer
}

func usage() {
	fmt.Print(
		"vcfFilter - Filter vcf records.\n\n" +
			"Usage:\n" +
			"  vcfFilter [options] input.vcf output.vcf\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

// criteria by which vcf records are filtered.
type criteria struct {
	chrom                          string
	groupFile                      string
	minPos                         int
	maxPos                         int
	minQual                        float64
	ref                            string
	alt                            []string
	biAllelicOnly                  bool
	substitutionsOnly              bool
	segregatingSitesOnly           bool
	removeNoAncestor               bool
	onlyPolarizableAncestors       bool
	weakToStrongOrStrongToWeakOnly bool
	noWeakToStrongOrStrongToWeak   bool
	refWeakAltStrongOnly           bool
	refStrongAltWeakOnly           bool
	notRefWeakAltStrong            bool
	notRefStrongAltWeak            bool
	id                             string
	formatExp                      string
	infoExp                        string
	includeMissingInfo             bool
	subSet                         float64
	minDaf                         float64
	maxDaf                         float64
	minDist                        int
}

// testingFuncs are a set of functions that must all return true to escape filter.
type testingFuncs []func(vcf.Vcf) bool

// passesTests runs all testingFuncs on a vcf record and returns true if all tests pass.
func passesTests(v vcf.Vcf, t testingFuncs) bool {
	for i := range t {
		if !t[i](v) {
			return false
		}
	}
	return true
}

// getTests parses the criteria struct to determine the testingFuncs.
func getTests(c criteria, header vcf.Header) testingFuncs {
	var answer testingFuncs

	if c.formatExp != "" {
		answer = append(answer, parseExpression(c.formatExp, header, true, c.includeMissingInfo)...) //raven's note: when tried to go run, got parseExpression undefined error. parseExpression is defined in cmd/vcfFilter/expression.go
	}

	if c.infoExp != "" {
		answer = append(answer, parseExpression(c.infoExp, header, false, c.includeMissingInfo)...)
	}

	if c.chrom != "" {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Chr == c.chrom
			})
	}

	if c.minPos != 0 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Pos >= c.minPos
			})
	}

	if c.maxPos != math.MaxInt64 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Pos <= c.maxPos
			})
	}

	if c.minDaf != 0 {
		if c.minDaf < 0 || c.minDaf > 1 {
			log.Fatalf("minDaf must be between 0 and 1.")
		}
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return popgen.VcfSampleDerivedAlleleFrequency(v) > c.minDaf
			})
	}

	if c.maxDaf != 1 {
		if c.maxDaf < 0 || c.maxDaf > 1 {
			log.Fatalf("maxDaf must be between 0 and 1.")
		}
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return popgen.VcfSampleDerivedAlleleFrequency(v) < c.maxDaf
			})
	}

	if c.maxDaf < c.minDaf {
		log.Fatalf("maxDaf must be less than minDaf.")
	}

	if c.minQual != 0 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Qual >= c.minQual
			})
	}

	if c.ref != "" {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Ref == c.ref
			})
	}

	if c.alt != nil {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if len(v.Alt) != len(c.alt) {
					return false
				}
				for i := range v.Alt {
					if v.Alt[i] != c.alt[i] {
						return false
					}
				}
				return true
			})
	}

	if c.biAllelicOnly {
		answer = append(answer, vcf.IsBiallelic)
	}

	if c.substitutionsOnly {
		answer = append(answer, vcf.IsSubstitution)
	}

	if c.segregatingSitesOnly {
		answer = append(answer, vcf.IsSegregating)
	}

	if c.removeNoAncestor {
		answer = append(answer, vcf.HasAncestor)
	}

	if c.onlyPolarizableAncestors {
		answer = append(answer, vcf.IsPolarizable)
	}
	if c.noWeakToStrongOrStrongToWeak {
		answer = append(answer, vcf.IsNotWeakToStrongOrStrongToWeak)
	}
	if c.weakToStrongOrStrongToWeakOnly {
		answer = append(answer, vcf.IsWeakToStrongOrStrongToWeak)
	}
	if c.refWeakAltStrongOnly {
		answer = append(answer, vcf.IsRefWeakAltStrong)
	}
	if c.refStrongAltWeakOnly {
		answer = append(answer, vcf.IsRefStrongAltWeak)
	}
	if c.notRefWeakAltStrong {
		answer = append(answer, vcf.IsNotRefWeakAltStrong)
	}
	if c.notRefStrongAltWeak {
		answer = append(answer, vcf.IsNotRefStrongAltWeak)
	}
	if c.id != "" {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				return v.Id == c.id
			})
	}
	if c.subSet < 1 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				r := rand.Float64()
				return r <= c.subSet
			})
	}
	return answer
}

func main() {
	var expectedNumArgs int = 2
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var groupFile *string = flag.String("groupFile", "", "Retains alleles from individuals in the input group file.")
	var minPos *int = flag.Int("minPos", 0, "Specifies the minimum position of the variant.")
	var maxPos *int = flag.Int("maxPos", math.MaxInt64, "Specifies the maximum position of the variant.")
	var minQual *float64 = flag.Float64("minQual", 0, "Specifies the minimum quality score.")
	var ref *string = flag.String("ref", "", "Specifies the reference field.")
	var alt *string = flag.String("alt", "", "Specifies the alt field.")
	var biAllelicOnly *bool = flag.Bool("biAllelicOnly", false, "Retains only biallelic variants in the output file. Not compatible with stdin.")
	var substitutionsOnly *bool = flag.Bool("substitutionsOnly", false, "Retains only substitution variants in the output file (removes INDEL variants).")
	var segregatingSitesOnly *bool = flag.Bool("segregatingSitesOnly", false, "Retains only variants that are segregating in at least one sample.")
	var removeNoAncestor *bool = flag.Bool("removeNoAncestor", false, "Retains only variants with an ancestor allele annotated in the info column.")
	var onlyPolarizableAncestors *bool = flag.Bool("onlyPolarizableAncestors", false, "Retains only variants that can be used to construct a derived allele frequency spectrum. Must have a subsitution where the ancestral allele matches either alt or ref.")
	var weakToStrongOrStrongToWeakOnly *bool = flag.Bool("weakToStrongOrStrongToWeakOnly", false, "Retains only variants that are weak to strong or strong to weak mutations.")
	var noWeakToStrongOrStrongToWeak *bool = flag.Bool("noWeakToStrongOrStrongToWeak", false, "Removes weak to strong variants and strong to weak variants.")
	var refWeakAltStrongOnly *bool = flag.Bool("refWeakAltStrongOnly", false, "Retains only variants that have a weak Ref allele and a strong Alt allele.")
	var refStrongAltWeakOnly *bool = flag.Bool("refStrongAltWeakOnly", false, "Retains only variants that have a strong Ref allele and a weak Alt allele.")
	var NotRefStrongAltWeak *bool = flag.Bool("notRefStrongAltWeak", false, "Removes variants that have a strong Ref alleles AND weak Alt alleles.")
	var NotRefWeakAltStrong *bool = flag.Bool("notRefWeakAltStrong", false, "Removes variants that have weak Ref allele AND a strong Alt allele.")
	var id *string = flag.String("id", "", "Specifies the rsID") //raven's note: added id string
	var formatExp *string = flag.String("format", "", "A logical expression (or a series of semicolon ';' delimited expressions) consisting of a tag and value present in the format field. Must be in double quotes (\"). "+
		"Expression can use the operators '>' '<' '=' '!=' '<=' '>'. For example, you can filter for variants with read depth greater than 100 and mapping quality greater or equal to 20 with the expression: \"DP > 100 ; MQ > 20\". "+
		"This tag is currently not supported for tags that have multiple values. When testing a vcf with multiple samples, the expression will only be tested on the first sample.")
	var infoExp *string = flag.String("info", "", "Identical to the 'format' tag, but tests the info field. The values of type 'Flag' in the info field"+
		"can be tested by including just the flag ID in the expression. E.g. To select all records with the flag 'GG' you would use the expression \"GG\".")
	var includeMissingInfo *bool = flag.Bool("includeMissingInfo", false, "When querying the records using the \"-info\" tag, include records where the queried tags are not present.")
	var subSet *float64 = flag.Float64("subSet", 1, "Proportion of variants to retain in output. Value must be between 0 and 1.")
	var minDaf *float64 = flag.Float64("minDaf", 0, "Set the minimum derived allele frequency for retained variants. Ancestral allele must be defined in INFO.")
	var maxDaf *float64 = flag.Float64("maxDaf", 1, "Set the maximum derived allele frequency for retained variants. Ancestral allele must be defined in INFO.")
	var minDist *int = flag.Int("minDistance", 0, "Remove variants that are within minDistance of another variant. File must be sorted by position.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	var altSlice []string
	if *alt != "" {
		altSlice = strings.Split(*alt, ",")
	}

	c := criteria{
		chrom:                          *chrom,
		minPos:                         *minPos,
		maxPos:                         *maxPos,
		minQual:                        *minQual,
		ref:                            *ref,
		alt:                            altSlice,
		biAllelicOnly:                  *biAllelicOnly,
		substitutionsOnly:              *substitutionsOnly,
		segregatingSitesOnly:           *segregatingSitesOnly,
		removeNoAncestor:               *removeNoAncestor,
		onlyPolarizableAncestors:       *onlyPolarizableAncestors,
		weakToStrongOrStrongToWeakOnly: *weakToStrongOrStrongToWeakOnly,
		noWeakToStrongOrStrongToWeak:   *noWeakToStrongOrStrongToWeak,
		refWeakAltStrongOnly:           *refWeakAltStrongOnly,
		refStrongAltWeakOnly:           *refStrongAltWeakOnly,
		notRefStrongAltWeak:            *NotRefStrongAltWeak,
		notRefWeakAltStrong:            *NotRefWeakAltStrong,
		id:                             *id,
		formatExp:                      *formatExp,
		infoExp:                        *infoExp,
		includeMissingInfo:             *includeMissingInfo,
		subSet:                         *subSet,
		minDaf:                         *minDaf,
		maxDaf:                         *maxDaf,
		minDist:                        *minDist,
	}

	var parseFormat, parseInfo bool
	if *formatExp != "" {
		parseFormat = true
	}
	if *infoExp != "" {
		parseInfo = true
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outfile := flag.Arg(1)

	total, removed := vcfFilter(infile, outfile, c, *groupFile, parseFormat, parseInfo, *setSeed)
	log.Printf("Processed  %d variants\n", total)
	log.Printf("Removed    %d variants\n", removed)
}
