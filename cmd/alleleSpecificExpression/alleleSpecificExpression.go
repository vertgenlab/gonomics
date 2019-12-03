package main

import (
	"flag"
	"fmt"
	//"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"os"
	//"strings"
)

func usage() {
	fmt.Print(
		"Allele Specific expression by SNPs\n1) Takes sam file of aligned reads\n2) Name of output files\n3) Vcf file containing SNPs\n")
	flag.PrintDefaults()
}

func chromAndPosToNumber(chrom int64, start int64) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func AlleleExpression(samFilename string, fileName string, v []*vcf.Vcf) {
	var aln *sam.SamAln = nil
	var done bool = false
	//var err error
	//sam file to read
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	fresh, _ := os.Create(fileName + "freshAllele.sam")
	defer fresh.Close()

	marine, _ := os.Create(fileName + "marineAllele.sam")
	defer marine.Close()

	sam.WriteHeaderToFileHandle(fresh, header)
	sam.WriteHeaderToFileHandle(marine, header)

	//chroms := chromInfo.SliceToMap(header.Chroms)

	var aMap = AlleleMap(v)

	var j, fSnpCount, mSnpCount int
	var currRefPos, currQueryPos, k int64
	var ok bool
	//var numDiscardReads int64
	//var code uint64
	var loc Location
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {

		fSnpCount, mSnpCount = 0, 0
		//use 1-based
		currRefPos = aln.Pos-1
		currQueryPos = 0

		//faIdx = chroms[aln.RName].Order

		for j = 0; j < len(aln.Cigar); j++ {
			switch aln.Cigar[j].Op {
			case 'S':
				//currRefPos += aln.Cigar[j].RunLength
				currQueryPos += aln.Cigar[j].RunLength
			case 'I':
				currQueryPos += aln.Cigar[j].RunLength
			case 'D':
				currRefPos += aln.Cigar[j].RunLength
			case 'M':
				for k = 0; k < aln.Cigar[j].RunLength; k++ {
					//code = chromAndPosToNumber(faIdx, currRefPos+k)
					loc = Location{Chr: aln.RName, Pos: currRefPos+k}
					_, ok = aMap[loc]
					if ok {
						//fmt.Println("Ref: ", dna.BasesToString(fMap[faIdx|currRefPos+k]), "Alt: ", dna.BasesToString(mMap[faIdx|currRefPos+k]))
						if aln.Seq[currQueryPos+k] == aMap[loc].Ref {
							fSnpCount++
						}
						if aln.Seq[currQueryPos+k] == aMap[loc].Alt {
							mSnpCount++
						}
					}
				}
				currRefPos += aln.Cigar[j].RunLength
				currQueryPos += aln.Cigar[j].RunLength
			}
		}
		if fSnpCount > mSnpCount {
			//fresh = append(fresh, aln)
			sam.WriteAlnToFileHandle(fresh, aln)
		}
		if mSnpCount > fSnpCount {
			//marine = append(marine, aln)
			sam.WriteAlnToFileHandle(marine, aln)
		}
	}
}

func AlleleMap(v []*vcf.Vcf) map[Location]Alleles {
	//ref := make(map[int64][]dna.Base)
	//alt := make(map[int64][]dna.Base)
	aMap := make(map[Location]Alleles)
	//var curr *vcf.Vcf
	var refRune, altRune []rune
	//var faIdx int64
	var loc Location
	//var ok bool
	for i := 0; i < len(v); i++ {
		curr := v[i]
		//faIdx = chromSize[curr.Chr].Order
		//_, ok = aMap[faIdx|curr.Pos-1]
		refRune = []rune(curr.Ref)
		altRune = []rune(curr.Alt)
		if len(refRune) == 1 && len(altRune) == 1 {
			if dna.RuneToBase(refRune[0]) != dna.Dot && dna.RuneToBase(altRune[0]) != dna.Dot {
				//code = chromAndPosToNumber(faIdx, curr.Pos-1)
				loc = Location{Chr: curr.Chr, Pos: curr.Pos-1}
				aMap[loc] = Alleles{Ref: dna.RuneToBase(refRune[0]), Alt: dna.RuneToBase(altRune[0])}
			}
		} else {

		}
	}
	return aMap
}

type Alleles struct {
	Ref dna.Base
	Alt dna.Base
}

type Location struct {
	Chr string
	Pos int64
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		//log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		log.Fatalf("\n\n./alleleSpecificExpression $sam $name $vcf\n\n")
	}
	vcfFile := vcf.Read(flag.Arg(2))
	//vcfFile := vcf.Read("snps.vcf")

	AlleleExpression(flag.Arg(0), flag.Arg(1), vcfFile)
}
