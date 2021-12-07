// Command Group: "VCF Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfInfo(filename string, outFile string, printNumDivergent bool) {
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

	v, _ := vcf.GoReadToChan(filename)
	out := fileio.EasyCreate(outFile)

	for current := range v {
		if printNumDivergent {
			if !vcf.HasAncestor(current) {
				log.Fatalf("Divergence can only be evaluated for VCF files with annotated ancestral alleles.")
			}
			if vcf.IsAltAncestor(current) {
				numDivergent++
			} else {
				numNotDivergent++
			}
		}
		if current.Ref == "A" {
			switch current.Alt[0] {
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
		} else if current.Ref == "C" {
			switch current.Alt[0] {
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
		} else if current.Ref == "G" {
			switch current.Alt[0] {
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
		} else if current.Ref == "T" {
			switch current.Alt[0] {
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
		} else if current.Ref == "N" {
			switch current.Alt[0] {
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
		} else if current.Ref == "-" {
			switch current.Alt[0] {
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
	}

	_, err = fmt.Fprintf(out, "Variant statistics on file: %s\n\n", filename)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Transitions\nA to G: %d. G to A: %d. C to T: %d. T to C: %d.\n\n", AtoG, GtoA, CtoT, TtoC)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Transversions\nA to C: %d. C to A: %d. G to T: %d. T to G: %d. A to T: %d. T to A: %d. C to G: %d. G to C: %d.\n\n", AtoC, CtoA, GtoT, TtoG, AtoT, TtoA, CtoG, GtoC)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Gaps Introduced\nA to Gap: %d. G to Gap: %d. C to Gap: %d. T to Gap: %d. N to Gap: %d.\n\n", AtoGap, GtoGap, CtoGap, TtoGap, NtoGap)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Gaps resolved\nGap to A: %d. Gap to C: %d. Gap to T: %d. Gap To G: %d. Gap to N: %d.\n\n", GapToA, GapToC, GapToT, GapToG, GapToN)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "N's introduced\nA to N: %d. T to N: %d. G to N: %d. C to N: %d.\n\n", AtoN, TtoN, GtoN, CtoN)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "N's resolved\nN to A: %d. N to G: %d. N to T: %d. N to C: %d.\n", NtoA, NtoG, NtoT, NtoC)
	exception.PanicOnErr(err)

	if printNumDivergent {
		_, err = fmt.Fprintf(out, "Number of Divergent Sites: %v. Number of non-divergent sites: %v.\n", numDivergent, numNotDivergent)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"vcfInfo - Provides summary statistics on an input VCF file.\n" +
			"Usage:\n" +
			"vcfInfo file.vcf out.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	var printNumDivergent *bool = flag.Bool("printNumDivergent", false, "Parse the ancestral allele information and return the number of divergent and non-divergent sites in the file.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	outFile := flag.Arg(1)

	vcfInfo(infile, outFile, *printNumDivergent)
}
