// Command Group: "Variant Calling & Annotation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
	"log"
	"strings"
)

func usage() {
	fmt.Print(
		"vcfContext - Determine the nucleotide context of mutations. \n" +
			"Currently only works for single-nucleotide variants (e.g. T > C becomes ATG > ACG).\n" +
			"Usage:\n" +
			" vcfContext [options] -r ref.fa -i in.vcf > out.txt\n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var ref *string = flag.String("r", "", "Reference fasta file (must be indexed). Required for -pad > 0.")
	var input *string = flag.String("i", "", "Input vcf file.")
	var output *string = flag.String("o", "stdout", "Output tsv file.")
	var includeComplements *bool = flag.Bool("includeComplements", false, "Do not merge totals from complementary mutations (e.g. separate C > A and G > T counts). Should only be used when looking at single-stranded mutations.")
	var pad *int = flag.Int("pad", 1, "Number of context bases on either side of variant (e.g. 0 == T, 1 == ATG, 2 == TATGA, ...")
	var verbose *int = flag.Int("v", 0, "Verbose output by setting to >0.")
	flag.Parse()
	flag.Usage = usage

	if *input == "" {
		usage()
		log.Fatalf("ERROR: must input a vcf")
	}

	if *ref == "" && *pad > 0 {
		usage()
		log.Fatalf("ERROR: must input a reference fasta for -pad > 0.")
	}

	vcfContext(*ref, *input, *output, *pad, *includeComplements, *verbose)
}

func vcfContext(fastaFile, vcfFile, output string, pad int, includeComplements bool, verbose int) {
	var err error
	ref := fasta.NewSeeker(fastaFile, "")
	out := fileio.EasyCreate(output)

	// MAP STRUCTURE:
	// Top level key of map is mutation w/o context (e.g. "T>C")
	// Second level key is reference sequence with context (e.g. "ATG")
	// Second level value is count observed.
	// e.g. m["T>C"]["ATG"] == 1
	m := initMap(pad)

	vcfChan, _ := vcf.GoReadToChan(vcfFile)
	var topKey, botKey string
	var keyFound bool
	var seq []dna.Base
	var numIncluded int
	for v := range vcfChan {
		// exclude multiallelic and indels
		if !vcf.IsBiallelic(v) || !vcf.IsSubstitution(v) || v.Pos == 1 {
			if verbose > 2 { // skip indels
				log.Printf("skipping\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
			}
			continue
		}

		topKey = v.Ref + ">" + v.Alt[0]
		// exclude if invalid bases in key
		if _, keyFound = m[topKey]; !keyFound {
			if verbose > 1 {
				log.Printf("skipping invalid bases\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
			}
			continue
		}

		seq, err = fasta.SeekByName(ref, v.Chr, (v.Pos-1)-pad, (v.Pos-1)+pad+1)
		// exclude if ref sequence does not match
		if err != nil || seq[pad] != dna.StringToBase(v.Ref) {
			if verbose > 1 {
				log.Printf("skipping error fetching seq\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
			}
			continue
		}

		botKey = dna.BasesToString(seq)
		// exclude if invalid bases in key
		if _, keyFound = m[topKey][botKey]; !keyFound {
			if verbose > 1 {
				log.Printf("skipping invalid bases in context\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
			}
			continue
		}

		numIncluded++
		m[topKey][botKey]++
	}

	if !includeComplements {
		mergeComplements(m)
	}

	_, err = fmt.Fprintf(out, getOutput(m))
	exception.PanicOnErr(err)
	if verbose > 0 {
		log.Printf("Processed %d variants\n", numIncluded)
	}

	err = ref.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

// initialize map with all valid permutations of nucleotides for given pad size
func initMap(pad int) map[string]map[string]int {
	m := make(map[string]map[string]int)

	// fill in all top level keys
	b := "ACGT"
	var i, j int
	for i = 0; i < len(b); i++ {
		for j = 0; j < len(b); j++ {
			if i == j {
				continue
			}
			m[string(b[i])+">"+string(b[j])] = make(map[string]int)
		}
	}

	// fill in 2nd level keys
	var pfs []string // pfs == possible flanking sequences
	var fullContextSeq string
	permute(b, "", pad*2, &pfs) // generate all possible flanking sequences
	for key, val := range m {   // loop over keys added in previous loop
		for i := range pfs {
			fullContextSeq = pfs[i][:pad] + string(key[0]) + pfs[i][pad:] // sandwich reference base (key[0]) between possible flanking sequence to generate all possible seqs for given ref basek
			val[fullContextSeq] = 0                                       // enter context sequence in map
		}
	}
	return m
}

// permute generates all possible permutations of characters present in b of length k and stores them in ans
func permute(b string, s string, k int, ans *[]string) {
	if k == 0 {
		*ans = append(*ans, s)
		return
	}

	for i := 0; i < len(b); i++ {
		permute(b, s+string(b[i]), k-1, ans)
	}
}

func mergeComplements(m map[string]map[string]int) {
	merge(m["C>A"], m["G>T"])
	delete(m, "G>T")
	merge(m["C>G"], m["G>C"])
	delete(m, "G>C")
	merge(m["C>T"], m["G>A"])
	delete(m, "G>A")
	merge(m["T>A"], m["A>T"])
	delete(m, "A>T")
	merge(m["T>C"], m["A>G"])
	delete(m, "A>G")
	merge(m["T>G"], m["A>C"])
	delete(m, "A>C")
}

func merge(m1, m2 map[string]int) {
	var b []dna.Base
	var keyFound bool
	for key := range m1 {
		b = dna.StringToBases(key)
		dna.ReverseComplement(b)
		if _, keyFound = m2[dna.BasesToString(b)]; !keyFound {
			log.Panicf("something went horribly wrong\n could not find key %s in map below\n%v\n", dna.BasesToString(b), m2)
		}
		m1[key] += m2[dna.BasesToString(b)]
	}
}

func getOutput(m map[string]map[string]int) string {
	var lines []string
	var k1, k2 string
	var val int
	var subMap map[string]int
	for k1, subMap = range m {
		for k2, val = range subMap {
			lines = append(lines, fmt.Sprintf("%s\t%s\t%d", k1, k2, val))
		}
	}
	slices.Sort(lines)
	return "Variant\tContext\tCount\n" + strings.Join(lines, "\n") + "\n"
}
