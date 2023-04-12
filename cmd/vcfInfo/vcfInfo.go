// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"golang.org/x/exp/slices"
)

func usage() {
	fmt.Print(
		"vcfInfo - Provides summary statistics on an input VCF file.\n" +
			"Usage:\n" +
			"vcfInfo -i file.vcf -types output1.txt -divergence output2.txt -context output3.txt -r reference.fasta\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var input *string = flag.String("i", "", "Input VCF file.")
	var types *string = flag.String("types", "", "Output file for summary info on the types of variants present in input file.")
	var divergence *string = flag.String("divergence", "", "Output file for ancestral divergence info. Parses the ancestral allele information from VCF record and returns the number of divergent and non-divergent sites.")
	var context *string = flag.String("context", "", "Output file for the nucleotide context of mutations. Currently only works for single-nucleotide variants (e.g. T > C becomes ATG > ACG).")
	var includeComplements *bool = flag.Bool("includeComplements", false, "For context mode only. Do not merge totals from complementary mutations (e.g. separate C > A and G > T counts). Should only be used when looking at single-stranded mutations.")
	var pad *int = flag.Int("pad", 1, "For context mode only. Number of context bases on either side of variant (e.g. 0 == T, 1 == ATG, 2 == TATGA, ...")
	var verbose *int = flag.Int("v", 0, "Verbose output by setting to >0.")
	var ref *string = flag.String("r", "", "For context mode only. Reference fasta file (must be indexed). Required for -pad > 0.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *input == "" {
		usage()
		log.Fatalln("ERROR: must define input vcf file with -i")
	}

	if *types == "" && *divergence == "" && *context == "" {
		usage()
		log.Fatalln("ERROR: must use at least one of -types, -divergence, or -context")
	}

	if *context != "" && *ref == "" && *pad > 0 {
		usage()
		log.Fatalf("ERROR: must input a reference fasta for context mode with -pad > 0.")
	}

	handleInputs(*input, *types, *divergence, *context, *includeComplements, *pad, *ref, *verbose)
}

func handleInputs(input, types, divergence, context string, includeComplements bool, pad int, fastaFile string, verbose int) {
	var typesOut, divergenceOut, contextOut *fileio.EasyWriter
	var ref *fasta.Seeker
	if types != "" {
		typesOut = fileio.EasyCreate(types)
		defer cleanup(typesOut)
	}
	if divergence != "" {
		divergenceOut = fileio.EasyCreate(divergence)
		defer cleanup(divergenceOut)
	}
	if context != "" {
		contextOut = fileio.EasyCreate(context)
		defer cleanup(contextOut)
		if pad > 0 {
			ref = fasta.NewSeeker(fastaFile, "")
			defer cleanup(ref)
		}
	}

	inputChan, _ := vcf.GoReadToChan(input)
	vcfInfo(input, inputChan, typesOut, divergenceOut, contextOut, ref, includeComplements, pad, verbose)
}

func vcfInfo(inputFile string, inputChan <-chan vcf.Vcf, typesOut, divergenceOut, contextOut *fileio.EasyWriter, ref *fasta.Seeker, includeComplements bool, pad, verbose int) {
	var AtoN, TtoN, GtoN, CtoN int
	var NtoA, NtoT, NtoG, NtoC int
	var AtoG, AtoT, AtoC, AtoGap int
	var TtoA, TtoC, TtoG, TtoGap int
	var GtoA, GtoC, GtoT, GtoGap int
	var CtoA, CtoT, CtoG, CtoGap int
	var GapToA, GapToC, GapToT, GapToG int
	var NtoGap, GapToN int
	var numDivergent, numNotDivergent int = 0, 0
	var err error

	var m map[string]map[string]int
	if contextOut != nil {
		// MAP STRUCTURE:
		// Top level key of map is mutation w/o context (e.g. "T>C")
		// Second level key is reference sequence with context (e.g. "ATG")
		// Second level value is count observed.
		// e.g. m["T>C"]["ATG"] == 1
		m = initMap(pad)
	}

	for v := range inputChan {
		switch v.Ref {
		case "A":
			switch v.Alt[0] {
			case "N":
				AtoN++
			case "T":
				AtoT++
			case "G":
				AtoG++
			case "C":
				AtoC++
			case "-":
				AtoGap++
			}
		case "C":
			switch v.Alt[0] {
			case "N":
				CtoN++
			case "T":
				CtoT++
			case "G":
				CtoG++
			case "A":
				CtoA++
			case "-":
				CtoGap++
			}
		case "G":
			switch v.Alt[0] {
			case "N":
				GtoN++
			case "C":
				GtoC++
			case "T":
				GtoT++
			case "A":
				GtoA++
			case "-":
				GtoGap++
			}
		case "T":
			switch v.Alt[0] {
			case "N":
				TtoN++
			case "A":
				TtoA++
			case "G":
				TtoG++
			case "C":
				TtoC++
			case "-":
				TtoGap++
			}
		case "N":
			switch v.Alt[0] {
			case "T":
				NtoT++
			case "G":
				NtoG++
			case "C":
				NtoC++
			case "A":
				NtoA++
			case "-":
				NtoGap++
			}
		case "-":
			switch v.Alt[0] {
			case "A":
				GapToA++
			case "T":
				GapToT++
			case "G":
				GapToG++
			case "C":
				GapToC++
			case "N":
				GapToN++
			}
		}

		if divergenceOut != nil {
			if !vcf.HasAncestor(v) {
				log.Fatalf("Divergence can only be evaluated for VCF files with annotated ancestral alleles.")
			}
			if vcf.IsAltAncestor(v) {
				numDivergent++
			} else {
				numNotDivergent++
			}
		}

		if contextOut != nil {
			vcfContext(v, m, ref, pad, verbose)
		}
	}

	if typesOut != nil {
		_, err = fmt.Fprintf(typesOut, "Variant statistics on file:\t%s\n\n", inputFile)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "Transitions\nA to G:\t%d\nG to A:\t%d\nC to T:\t%d\nT to C:\t%d\n\n", AtoG, GtoA, CtoT, TtoC)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "Transversions\nA to C:\t%d\nC to A:\t%d\nG to T:\t%d\nT to G:\t%d\nA to T:\t%d\nT to A:\t%d\nC to G:\t%d\nG to C:\t%d\n\n", AtoC, CtoA, GtoT, TtoG, AtoT, TtoA, CtoG, GtoC)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "Gaps Introduced\nA to Gap:\t%d\nG to Gap:\t%d\nC to Gap:\t%d\nT to Gap:\t%d\nN to Gap:\t%d\n\n", AtoGap, GtoGap, CtoGap, TtoGap, NtoGap)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "Gaps resolved\nGap to A:\t%d\nGap to C:\t%d\nGap to T:\t%d\nGap To G:\t%d\nGap to N:\t%d\n\n", GapToA, GapToC, GapToT, GapToG, GapToN)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "N's introduced\nA to N:\t%d\nT to N:\t%d\nG to N:\t%d\nC to N:\t%d\n\n", AtoN, TtoN, GtoN, CtoN)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(typesOut, "N's resolved\nN to A:\t%d\nN to G:\t%d\nN to T:\t%d\nN to C:\t%d\n\n", NtoA, NtoG, NtoT, NtoC)
		exception.PanicOnErr(err)
	}

	if divergenceOut != nil {
		_, err = fmt.Fprintf(divergenceOut, "Variant statistics on file:\t%s\n\n", inputFile)
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(divergenceOut, "Number of Divergent Sites:\t%v\nNumber of non-divergent sites:\t%v\n", numDivergent, numNotDivergent)
		exception.PanicOnErr(err)
	}

	if contextOut != nil {
		if !includeComplements {
			mergeComplements(m)
		}
		_, err = fmt.Fprintf(contextOut, getOutput(m))
		exception.PanicOnErr(err)
	}
}

func cleanup(f io.Closer) {
	var err error
	err = f.Close()
	exception.PanicOnErr(err)
}

func vcfContext(v vcf.Vcf, m map[string]map[string]int, ref *fasta.Seeker, pad int, verbose int) {
	var topKey, botKey string
	var keyFound bool
	var seq []dna.Base
	var numIncluded int
	var err error

	// exclude multiallelic and indels
	if !vcf.IsBiallelic(v) || !vcf.IsSubstitution(v) || v.Pos == 1 {
		if verbose > 2 { // skip indels
			log.Printf("skipping\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	topKey = v.Ref + ">" + v.Alt[0]
	// exclude if invalid bases in key
	if _, keyFound = m[topKey]; !keyFound {
		if verbose > 1 {
			log.Printf("skipping invalid bases\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	if pad > 0 {
		seq, err = fasta.SeekByName(ref, v.Chr, (v.Pos-1)-pad, (v.Pos-1)+pad+1)
	} else {
		seq = dna.StringToBases(v.Ref)
	}
	// exclude if ref sequence does not match
	if err != nil || seq[pad] != dna.StringToBase(v.Ref) {
		if verbose > 1 {
			log.Printf("skipping error fetching seq\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	botKey = dna.BasesToString(seq)
	// exclude if invalid bases in key
	if _, keyFound = m[topKey][botKey]; !keyFound {
		if verbose > 1 {
			log.Printf("skipping invalid bases in context\t%s:%d\t%s>%s\n", v.Chr, v.Pos, v.Ref, strings.Join(v.Alt, ","))
		}
		return
	}

	numIncluded++
	m[topKey][botKey]++
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
