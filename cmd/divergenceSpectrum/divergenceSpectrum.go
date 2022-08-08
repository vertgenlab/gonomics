// Command Group: "Sequence Evolution & Reconstruction"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func divergenceSpectrum(inBed string, inVcf string, outFile string) {
	var overlappingIntervals []interval.Interval
	var currSpectrum []int
	var err error

	variants, _ := vcf.Read(inVcf)
	varIntervals := make([]interval.Interval, 0)
	for i := range variants {
		varIntervals = append(varIntervals, variants[i])
	}

	tree := interval.BuildTree(varIntervals)
	bedChan := bed.GoReadToChan(inBed)
	out := fileio.EasyCreate(outFile)

	for i := range bedChan {
		overlappingIntervals = interval.Query(tree, i, "any")
		currSpectrum = make([]int, 6)
		for j := range overlappingIntervals {
			currSpectrum[mutationType(overlappingIntervals[j].(vcf.Vcf))]++
		}
		i.Annotation = make([]string, 6) //replace whatever's in the input bed annotation with our output
		for k := range currSpectrum {
			i.Annotation[k] = fmt.Sprint(currSpectrum[k])
		}
		i.FieldsInitialized = 13
		bed.WriteBed(out, i)
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// mutationType returns an int corresponding to six categories of mutations. (0: A->G/T->C) (1: G->A/C->T) (2: A->T/T->A) (3: G->C/C->G) (4: A->C/T->G) (5: C->A/G->T)
// 0 and 1 are transitions, 2-5 are transversions.
func mutationType(v vcf.Vcf) int {
	switch v.Ref {
	case "A":
		switch v.Alt[0] {
		case "C":
			return 4
		case "G":
			return 0
		case "T":
			return 2
		default:
			log.Fatalf("Poorly formed VCF. Ref: %s. Alt: %s.", v.Ref, v.Alt[0])
		}
	case "C":
		switch v.Alt[0] {
		case "A":
			return 5
		case "G":
			return 3
		case "T":
			return 1
		default:
			log.Fatalf("Poorly formed VCF. Ref: %s. Alt: %s.", v.Ref, v.Alt[0])
		}
	case "G":
		switch v.Alt[0] {
		case "A":
			return 1
		case "C":
			return 3
		case "T":
			return 5
		default:
			log.Fatalf("Poorly formed VCF. Ref: %s. Alt: %s.", v.Ref, v.Alt[0])
		}
	case "T":
		switch v.Alt[0] {
		case "A":
			return 2
		case "C":
			return 0
		case "G":
			return 4
		default:
			log.Fatalf("Poorly formed VCF. Ref: %s. Alt: %s.", v.Ref, v.Alt[0])
		}
	default:
		log.Fatalf("Poorly formed VCF. Ref: %s. Alt: %s.", v.Ref, v.Alt[0])
	}

	return -1
}

func usage() {
	fmt.Print("divergenceSpectrum - Determine the mutation spectrum for divergent sites in each region of an input bed file.\n" +
		"Usage:\n" +
		"divergenceSpectrum input.bed divergentSites.vcf output.txt\n\n" +
		"Requires an input vcf file of all divergent sites between two genomes, which can be generated with multiFaToVcf." +
		"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inBed := flag.Arg(0)
	inVcf := flag.Arg(1)
	outFile := flag.Arg(2)

	divergenceSpectrum(inBed, inVcf, outFile)
}
