package main

import (
	"flag"
	"fmt"
	"strings"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"

	"log"
	//"os"
)

func usage() {
	fmt.Print(
		"alleleSplit - a tool to separate allele specific heterozygous alignments into two split sams\n" +
			"\nusage:\n" +
			"\t./alleleSplit [options] het.sam SNPs.vcf.gz\n\n" +
			"required:\n" +
			"\t\t\tAll names specified in x,a,b must match columns of SNPs.vcf.gz exactly\n" +
			"\t--x\t\tlist of f1 hybrid samples\n" + //that map to columns in SNPs.vcf.gz separate multiple samples by commas\n" +
			"\t--a\t\tname of allele one of cross\n" + //must match name found in SNPs.vcf.gz\n" +
			"\t--b\t\tname of allele two of cross\n\n" + //must match name found in SNPs.vcf.gz"

			"options:\n" +
			"\t--name\t\tprefix string to name output sam files\n" +
			"\t--filter\t[AS/medium/strict/wrong]\n\t\t\tExtra filters to apply to input SNPs\n\n" +
			"\t\t\tAS: default when het input is 1-2 samples all must have het genotype in vcf\n" +
			"\t\t\tmedium: het input cohort >2, requires at least one het genotype to be considered\n" +
			"\t\t\tstrict: het input cohort >2 requires every sample to have hets in vcf record\n" +
			"\t\t\twrong: bad data to calibrate VQSR filter\n\n" +
			"./alleleSplit --x wgsAxB,atacAxB --a wgsAa --b wgsBb --name atacAxB het.sam SNPs.vcf.gz\n\n" +
			"./alleleSplit --filter snps.vcf samples.csv [AS/medium/strict/wrong] output\n\n")
}
func main() {
	var cross *string = flag.String("x", "", "string containing hets")
	var alleleOne *string = flag.String("a", "", "allele one")
	var alleleTwo *string = flag.String("b", "", "allele two")
	var name *string = flag.String("name", "f1Hybrd", "prefix to name output files")
	var filter *bool = flag.Bool("filter", false, "filter vcf")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *filter {
		var expectedNumArgs int = 4
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
		}
		file := flag.Arg(0)
		sampleSheet := flag.Arg(1)
		vcfFilter := flag.Arg(2)
		output := fileio.EasyCreate(flag.Arg(3))
		header := vcf.GetHeader(file)

		//output := fileio.MustCreate(flag.Arg(3))

		defer output.Close()
		vcf.WriteHeaderDev(output, header)
		dict := vcf.HeaderToMaps(file)

		vcfPipe := make(chan *vcf.Vcf)

		AaAa, AAaa := vcf.ReadFilterList(sampleSheet, dict.HapIdx)
		index1, index2 := vcf.MapNameToIndex(dict.HapIdx, AaAa), vcf.MapNameToIndex(dict.HapIdx, AAaa)
		go vcf.ReadToChan(file, vcfPipe)
		logHeader := fmt.Sprintf("Chrom\tPos\tRef\tAlt\t%s\t%s", strings.Join(AaAa, "\t"), strings.Join(AAaa, "\t"))

		fmt.Printf("%s\n", logHeader)

		var group []vcf.Haplotype
		for gVcf := range vcfPipe {
			group = vcf.GenotypeHelper(gVcf)
			if vcf.ApplyFilter(vcfFilter, group, index1, index2) {
				vcf.WriteVcf(output, gVcf)
				vcf.PrettyShuffle(gVcf, index1, index2)
			}

		}
		//   for record := range vcfPipe {
		//    	if vcf.WrongFilter(record, vcf.MapNameToIndex(dict.HapIdx, AaAa), vcf.MapNameToIndex(dict.HapIdx, AAaa)) {
		//    		vcf.WriteVcf(output, record)
		//           	vcf.ViewGenotypeVcf(record)
		//    	}
		//    }

	} else {
		samfile := flag.Arg(0)
		snps := flag.Arg(1)
		catchError(*cross, *alleleOne, *alleleTwo)
		vcf.SnpSearch(samfile, snps, *cross, *alleleOne, *alleleTwo, *name)
	}
}

/*
file := flag.Arg(0)
	header := vcf.GetHeader(file)
	output := fileio.MustCreate()
	defer output.Close()
	vcf.WriteHeaderDev(output, header)
	dict := vcf.HeaderToMaps(file)
	vcfPipe := make(chan *vcf.Vcf)*/

func catchError(cross string, one string, two string) {
	if len(flag.Args()) < 2 {
		flag.Usage()
	}
	if strings.Compare(cross, "") == 0 || strings.Compare(one, "") == 0 || strings.Compare(two, "") == 0 {
		flag.Usage()
	}
}

/*
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
}*/
