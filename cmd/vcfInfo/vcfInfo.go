package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func vcfInfo(filename string) {
	var AtoN, TtoN, GtoN, CtoN int32
	var NtoA, NtoT, NtoG, NtoC int32
	var AtoG, AtoT, AtoC, AtoGap int32
	var TtoA, TtoC, TtoG, TtoGap int32
	var GtoA, GtoC, GtoT, GtoGap int32
	var CtoA, CtoT, CtoG, CtoGap int32
	var GaptoA, GaptoC, GaptoT, GaptoG int32
	var NtoGap, GaptoN int32

	v := vcf.Read(filename)
	fmt.Printf("len: %d", len(v))

	for _, current := range v {
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
				GaptoA++
			case "T":
				GaptoT++
			case "G":
				GaptoG++
			case "C":
				GaptoC++
			case "N":
				GaptoN++
			}
		}
	}

	fmt.Printf("Variant statistics on file: %s\n\n", filename)
	fmt.Printf("Transitions\n")
	fmt.Printf("A to G: %d. G to A: %d. C to T: %d. T to C: %d.\n\n", AtoG, GtoA, CtoT, TtoC)
	fmt.Printf("Transversions\n")
	fmt.Printf("A to C: %d. C to A: %d. G to T: %d. T to G: %d. A to T: %d. T to A: %d. C to G: %d. G to C: %d.\n\n", AtoC, CtoA, GtoT, TtoG, AtoT, TtoA, CtoG, GtoC)
	fmt.Printf("Gaps Introduced\n")
	fmt.Printf("A to Gap: %d. G to Gap: %d. C to Gap: %d. T to Gap: %d. N to Gap: %d.\n\n", AtoGap, GtoGap, CtoGap, TtoGap, NtoGap)
	fmt.Printf("Gaps resolved\n")
	fmt.Printf("Gap to A: %d. Gap to C: %d. Gap to T: %d. Gap To G: %d. Gap to N: %d.\n\n", GaptoA, GaptoC, GaptoT, GaptoG, GaptoN)
	fmt.Printf("N's introduced\n")
	fmt.Printf("A to N: %d. T to N: %d. G to N: %d. C to N: %d.\n\n", AtoN, TtoN, GtoN, CtoN)
	fmt.Printf("N's resolved\n")
	fmt.Printf("N to A: %d. N to G: %d. N to T: %d. N to C: %d.\n", NtoA, NtoG, NtoT, NtoC)
}

func usage() {
	fmt.Print(
		"vcfInfo - Provides summary statistics on an input VCF file.\n" +
			"Usage:\n" +
			"vcfInfo file.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)

	vcfInfo(infile)
}
