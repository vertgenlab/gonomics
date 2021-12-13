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

	_, err = fmt.Fprintf(out, "Variant statistics on file:\t%s\n\n", filename)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Transitions\nA to G:\t%d\nG to A:\t%d\nC to T:\t%d\nT to C:\t%d\n\n", AtoG, GtoA, CtoT, TtoC)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Transversions\nA to C:\t%d\nC to A:\t%d\nG to T:\t%d\nT to G:\t%d\nA to T:\t%d\nT to A:\t%d\nC to G:\t%d\nG to C:\t%d\n\n", AtoC, CtoA, GtoT, TtoG, AtoT, TtoA, CtoG, GtoC)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Gaps Introduced\nA to Gap:\t%d\nG to Gap:\t%d\nC to Gap:\t%d\nT to Gap:\t%d\nN to Gap:\t%d\n\n", AtoGap, GtoGap, CtoGap, TtoGap, NtoGap)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "Gaps resolved\nGap to A:\t%d\nGap to C:\t%d\nGap to T:\t%d\nGap To G:\t%d\nGap to N:\t%d\n\n", GapToA, GapToC, GapToT, GapToG, GapToN)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "N's introduced\nA to N:\t%d\nT to N:\t%d\nG to N:\t%d\nC to N:\t%d\n\n", AtoN, TtoN, GtoN, CtoN)
	exception.PanicOnErr(err)
	_, err = fmt.Fprintf(out, "N's resolved\nN to A:\t%d\nN to G:\t%d\nN to T:\t%d\nN to C:\t%d\n\n", NtoA, NtoG, NtoT, NtoC)
	exception.PanicOnErr(err)

	if printNumDivergent {
		_, err = fmt.Fprintf(out, "Number of Divergent Sites:\t%v\nNumber of non-divergent sites:\t%v\n", numDivergent, numNotDivergent)
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
