package main

import (
	"flag"
	"fmt"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	//"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//"sync"
)

func usage() {
	fmt.Print(
		"alleleSplit - a tool that separates a heterozygous sam alignment into different alignments by alleles\n" +
			"Usage:\n" +
			"./alleleSplit [options] input.sam input.vcf\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	var parental arrayFlags
	flag.Var(&parental, "parent", "Define of two parential that appears homozygous in the genotype Vcf. Include second `name -parent two -f1 name`")
	var f1Genome *string = flag.String("f1", "", "F1 hybrid sample that appears heterozygous in genotype Vcf. Include `name -parent one -parent two`")
	var sampleName *bool = flag.Bool("samples", false, "Get names of samples that appear in Vcf header. (Default: `/dev/stdout`)")

	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("\n\nError: unexpected number of arguments...\n\n")
	}
	if *sampleName {
		file := fileio.EasyOpen(flag.Arg(0))
		defer file.Close()
		header := vcf.ReadHeader(file)
		fmt.Printf("%s", vcf.PrintSampleNames(header))
	} else {
		if len(parental) != 2 {
			log.Fatalf("Error: must provide exactly 2 parental genomes\n")
		}
		SnpSearch(flag.Arg(0), flag.Arg(1), *f1Genome, parental[0], parental[1], *f1Genome)
	}

	//BasicAlleleExpression(flag.Arg(0), vcf.Read(flag.Arg(1)), flag.Arg(2), *refPrefix, *altPrefix)
}

//Define flag value as an array. Used to define parent genomes.
type arrayFlags []string

func (i *arrayFlags) String() string {
	return "my string representation"
}
func (i *arrayFlags) Set(value string) error {
	*i = append(*i, value)
	return nil
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

func BasicAlleleExpression(samFilename string, v []*vcf.Vcf, fileName string, ref string, alt string) {
	var aln *sam.SamAln = nil
	var done bool = false
	//var err error
	//sam file to read
	samFile := fileio.EasyOpen(samFilename)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	refFile := fileio.EasyCreate(fileName + "_" + ref + "Allele.sam")
	defer refFile.Close()

	altFile := fileio.EasyCreate(fileName + "_" + alt + "Allele.sam")
	defer altFile.Close()
	sam.WriteHeaderToFileHandle(refFile, header)
	sam.WriteHeaderToFileHandle(altFile, header)
	//var un *fileio.EasyWriter = nil
	//if unFlag {
	//	un = fileio.EasyCreate(fileName + "_undetermined.sam")
	//	defer un.Close()
	//	sam.WriteHeaderToFileHandle(un, header)
	//}
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
			//} else if un != nil {
			//sam.WriteAlnToFileHandle(un, aln)
		} else {
			//Skip read
		}
	}
}*/

/*
func GoRoutinesSnpSearch(samfile string, genotypeVcf string, parentOne string, parentTwo, prefix string, f1 string, threads int) {
	var wg sync.WaitGroup
	gvcf := make(chan *vcf.Vcf)
	genotypeReader := fileio.EasyOpen(genotypeVcf)
	defer genotypeReader.Close()

	vcfHeader := vcf.ReadHeader(genotypeReader)

	sampleHash := vcf.HeaderToMaps(vcfHeader)

	go vcf.ReadToChan(genotypeReader, gvcf)

	//children := strings.Split(f1, ",")
	//parents := []string{parentOne, parentTwo}
	//hets, homs := vcf.MapNameToIndex(sampleHash.Index, children), vcf.MapNameToIndex(sampleHash.Index, parents)

	snpDb := make(map[uint64]*vcf.GVcf)
	var parentalOne, parentalTwo, fOne int16 = sampleHash.IndexAllele[parentOne], sampleHash.IndexAllele[parentTwo], sampleHash.IndexAllele[f1]
	for genotype := range gvcf {
		if vcf.ASFilter(genotype, parentalOne, parentalTwo, fOne) {
			snpDb = vcf.GenotypeToMap(genotype, sampleHash.FaIndex)
		}
	}
	samFile := fileio.EasyOpen(samfile)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	samReader := make(chan *sam.SamAln)
	go sam.ReadToChan(samFile, samReader)
	wg.Wait()
	var wgReader, wgWriter sync.WaitGroup
	childOne := make(chan *sam.SamAln)
	childTwo := make(chan *sam.SamAln)

	parents := []string{parentOne, parentTwo}
	for i := 0; i < threads; i++ {
		wgReader.Add(1)
		go snpAnalysis(snpDb, sampleHash, parents, samReader, childOne, childTwo, &wgReader)
	}
	wgWriter.Add(2)
	go sam.SamChanToFile(childOne, fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentOne), header, &wgWriter)
	go sam.SamChanToFile(childTwo, fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentTwo), header, &wgWriter)
	wgReader.Wait()
	close(childOne)
	close(childTwo)
	wgWriter.Wait()
}

func snpAnalysis(snpDb map[uint64]*vcf.GVcf, sampleHash *vcf.SampleIdMap, parents []string, samReader <-chan *sam.SamAln, childOne chan<- *sam.SamAln, childTwo chan<- *sam.SamAln, wg *sync.WaitGroup) {
	//for read, done := sam.NextAlignment(samFile); done != true; read, done = sam.NextAlignment(samFile) {
	for read := range samReader {
		parentAllele1, parentAllele2 := 0, 0
		var target int64 = read.Pos - 1
		var query int64 = 0
		var code uint64
		var ok bool
		var gV *vcf.GVcf
		for i := 0; i < len(read.Cigar); i++ {
			switch read.Cigar[i].Op {
			case 'S':
				query += read.Cigar[i].RunLength
			case 'I':
				//code = vcf.ChromPosToUInt64(int(sampleHash.FaIndex[read.RName]), int(target))
				//_, ok = snpDb[code]
				//if ok {
				//	gV = snpDb[code]
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleTwo]) == 0 {
				//		parentAllele1++
				//	}
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleTwo]) == 0 {
				//		parentAllele2++
				//	}
				//}
				query += read.Cigar[i].RunLength
			case 'D':
				code = vcf.ChromPosToUInt64(int(sampleHash.FaIndex[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				var j int64
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = vcf.ChromPosToUInt64(int(sampleHash.FaIndex[read.RName]), int(target+j))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleTwo]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleTwo]) == 0 {
							parentAllele2++
						}
					}

				}
				target += read.Cigar[i].RunLength
				query += read.Cigar[i].RunLength
			}
		}
		if parentAllele1 > parentAllele2 {
			childOne <- read
			//sam.WriteAlnToFileHandle(childOne, read)
		} else if parentAllele2 > parentAllele1 {
			//sam.WriteAlnToFileHandle(childTwo, read)
			childTwo <- read
		} else {
			//Skip read
		}
	}
	wg.Done()
}*/
