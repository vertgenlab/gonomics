package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
)

func usage() {
	fmt.Print(
		"alleleSplit - a tool that separates a heterozygous sam alignment into different alleles\n" +
			"1) Takes sam file of aligned reads\n2)Vcf file containing SNPs\n3)Name of output files\n\n" +
			"usage:\n" +
			"./alleleSplit [options] input.sam input.vcf output.sam\n" +
			"\t--ref ref.output\n" +
			"\t\tname output ref allele sam file\n" +
			"\t--alt alt.output\n" +
			"\t\tname output alt allele sam file\n" +
			"\t--undetermined undetermined.sam\n" +
			"\t\toutputs third file containing reads discarded from analysis\n")
}

type Alleles struct {
	Ref dna.Base
	Alt dna.Base
}

type Location struct {
	Chr string
	Pos int64
}

func AlleleMap(v []*vcf.Vcf) map[Location]Alleles {
	aMap := make(map[Location]Alleles)
	var refRune, altRune []rune
	var loc Location
	for i := 0; i < len(v); i++ {
		curr := v[i]
		refRune = []rune(curr.Ref)
		altRune = []rune(curr.Alt)
		if len(refRune) == 1 && len(altRune) == 1 {
			if dna.RuneToBase(refRune[0]) != dna.Dot && dna.RuneToBase(altRune[0]) != dna.Dot {
				loc = Location{Chr: curr.Chr, Pos: curr.Pos - 1}
				aMap[loc] = Alleles{Ref: dna.RuneToBase(refRune[0]), Alt: dna.RuneToBase(altRune[0])}
			}
		} else {

		}
	}
	return aMap
}

func AlleleExpression(samFilename string, v []*vcf.Vcf, fileName string, ref string, alt string, unFlag bool) {
	var aln *sam.SamAln = nil
	var done bool = false
	//var err error
	//sam file to read
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	refFile, _ := os.Create(fileName + "_" + ref + "Allele.sam")
	defer refFile.Close()

	altFile, _ := os.Create(fileName + "_" + alt + "Allele.sam")
	defer altFile.Close()
	sam.WriteHeaderToFileHandle(refFile, header)
	sam.WriteHeaderToFileHandle(altFile, header)
	var un *os.File = nil
	if unFlag {
		un, _ := os.Create(fileName + "_undetermined.sam")
		defer un.Close()
		sam.WriteHeaderToFileHandle(un, header)
	}
	var aMap = AlleleMap(v)

	var j, refCount, altCount int
	var currRefPos, currQueryPos, k int64
	var ok bool
	var loc Location
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		refCount, altCount = 0, 0
		currRefPos = aln.Pos - 1
		currQueryPos = 0
		for j = 0; j < len(aln.Cigar); j++ {
			switch aln.Cigar[j].Op {
			case 'S':
				currQueryPos += aln.Cigar[j].RunLength
			case 'I':
				currQueryPos += aln.Cigar[j].RunLength
			case 'D':
				currRefPos += aln.Cigar[j].RunLength
			case 'M':
				for k = 0; k < aln.Cigar[j].RunLength; k++ {
					loc = Location{Chr: aln.RName, Pos: currRefPos + k}
					_, ok = aMap[loc]
					if ok {
						if aln.Seq[currQueryPos+k] == aMap[loc].Ref {
							refCount++
						}
						if aln.Seq[currQueryPos+k] == aMap[loc].Alt {
							altCount++
						}
					}
				}
				currRefPos += aln.Cigar[j].RunLength
				currQueryPos += aln.Cigar[j].RunLength
			}
		}
		if refCount > altCount {
			sam.WriteAlnToFileHandle(refFile, aln)
		} else if altCount > refCount {
			sam.WriteAlnToFileHandle(altFile, aln)
		} else if un != nil {
			sam.WriteAlnToFileHandle(un, aln)
		} else {
			//Skip read
		}
	}
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)

	var refPrefix *string = flag.String("ref", "ref", "name of ref allele")
	var altPrefix *string = flag.String("alt", "alt", "name of alt allele")
	var undetermined *bool = flag.Bool("undetermined", false, "undetermined.sam")
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("\n\n./alleleSplit [options] $sam $vcf $name\n\n")
	}
	AlleleExpression(flag.Arg(0), vcf.Read(flag.Arg(1)), flag.Arg(2), *refPrefix, *altPrefix, *undetermined)
}
