// Command Group: "VCF Tools"

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

func vcfFilter(infile string, outfile string, c criteria, groupFile string, parseFormat bool, parseInfo bool) (total, removed int){
	records, header := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(outfile)
	tests := getTests(c, header)

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

	for v := range records {
		total++
		if groupFile != "" {
			v.Samples = filterRecordsSamplesToKeep(v.Samples, samplesToKeep)
		}

		if parseFormat {
			v = vcf.ParseFormat(v, header)
		}

		if parseInfo {
			v = vcf.ParseInfo(v, header)
		}

		if !passesTests(v, tests) {
			removed++
			continue
		}

		vcf.WriteVcf(out.File, v)
	}

	err := out.Close()
	if err != nil {
		log.Panic(err)
	}
	return
}

func filterRecordsSamplesToKeep(recordSamples []vcf.GenomeSample, samplesToKeep []int) []vcf.GenomeSample {
	var answer []vcf.GenomeSample
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
	chrom                    string
	groupFile                string
	minPos                   int
	maxPos                   int
	minQual                  float64
	ref                      string
	alt                      []string
	biAllelicOnly            bool
	substitutionsOnly        bool
	segregatingSitesOnly     bool
	removeNoAncestor         bool
	onlyPolarizableAncestors bool
	formatExp                string
	infoExp                  string
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
		answer = append(answer, parseExpression(c.formatExp, header, true)...)
	}

	if c.infoExp != "" {
		answer = append(answer, parseExpression(c.infoExp, header, false)...)
	}

	if c.chrom != "" {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if v.Chr != c.chrom {
					return false
				}
				return true
			})
	}

	if c.minPos != 0 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if v.Pos < c.minPos {
					return false
				}
				return true
			})
	}

	if c.maxPos != math.MaxInt64 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if v.Pos > c.maxPos {
					return false
				}
				return true
			})
	}

	if c.minQual != 0 {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if v.Qual < c.minQual {
					return false
				}
				return true
			})
	}

	if c.ref != "" {
		answer = append(answer,
			func(v vcf.Vcf) bool {
				if v.Ref != c.ref {
					return false
				}
				return true
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
	return answer
}

func main() {
	var expectedNumArgs int = 2
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var groupFile *string = flag.String("groupFile", "", "Retains alleles from individuals in the input group file.")
	var minPos *int = flag.Int("minPos", 0, "Specifies the minimum position of the variant.")
	var maxPos *int = flag.Int("maxPos", math.MaxInt64, "Specifies the maximum position of the variant.")
	var minQual *float64 = flag.Float64("minQual", 0, "Specifies the minimum quality score.")
	var ref *string = flag.String("ref", "", "Specifies the reference field.")
	var alt *string = flag.String("alt", "", "Specifies the alt field.")
	var biAllelicOnly *bool = flag.Bool("biAllelicOnly", false, "Retains only biallelic variants in the output file.")
	var substitutionsOnly *bool = flag.Bool("substitutionsOnly", false, "Retains only substitution variants in the output file (removes INDEL variants).")
	var segregatingSitesOnly *bool = flag.Bool("segregatingSitesOnly", false, "Retains only variants that are segregating in at least one sample.")
	var removeNoAncestor *bool = flag.Bool("removeNoAncestor", false, "Retains only variants with an ancestor allele annotated in the info column.")
	var onlyPolarizableAncestors *bool = flag.Bool("onlyPolarizableAncestors", false, "Retains only variants that can be used to construct a derived allele frequency spectrum. Must have a subsitution where the ancestral allele matches either alt or ref.")
	var formatExp *string = flag.String("format", "", "A logical expression (or a series of semicolon ';' delimited expressions) consisting of a tag and value present in the format field. Must be in double quotes (\"). "+
		"Expression can use the operators '>' '<' '=' '!=' '<=' '>'. For example, you can filter for variants with read depth greater than 100 and mapping quality greater or equal to 20 with the expression: \"DP > 100 ; MQ > 20\". "+
		"This tag is currently not supported for tags that have multiple values. When testing a vcf with multiple samples, the expression will only be tested on the first sample.")
	var infoExp *string = flag.String("info", "", "Identical to the 'format' tag, but tests the info field. The values of type 'Flag' in the info field"+
		"can be tested by including just the flag ID in the expression. E.g. To select all records with the flag 'GG' you would use the expression \"GG\".")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	var altSlice []string
	if *alt != "" {
		altSlice = strings.Split(*alt, ",")
	}

	c := criteria{
		chrom:                    *chrom,
		minPos:                   *minPos,
		maxPos:                   *maxPos,
		minQual:                  *minQual,
		ref:                      *ref,
		alt:                      altSlice,
		biAllelicOnly:            *biAllelicOnly,
		substitutionsOnly:        *substitutionsOnly,
		segregatingSitesOnly:     *segregatingSitesOnly,
		removeNoAncestor:         *removeNoAncestor,
		onlyPolarizableAncestors: *onlyPolarizableAncestors,
		formatExp:                *formatExp,
		infoExp:                  *infoExp,
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

	total, removed := vcfFilter(infile, outfile, c, *groupFile, parseFormat, parseInfo)
	fmt.Printf("Processed  %d variants\nRemoved    %d variants\n", total, removed)
}
