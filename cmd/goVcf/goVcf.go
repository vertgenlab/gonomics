package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"sync"
	//"strings"
)

func usage() {
	fmt.Print(
		"goVcf - a gonomics tool for analyzing genotyped vcfs\n\n" +
			"Usage:\n" +
			"  ./goVcf [options] .vcf/.gz\n\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	var filterAs *bool = flag.Bool("filter", false, "Filter for Vcf SNPs regions where the provided parental genomes are homozgous and F1 is heterozygous\nMust define two parental genomes and one F1 hybrid")
	var parental arrayFlags
	flag.Var(&parental, "parent", "Define of two parential that appears homozygous in the genotype Vcf\nMust match Vcf Header exactly``")
	var f1Genome *string = flag.String("f1", "", "Define of f1 hybrid sample that appears heterozygous in the genotype Vcf,\nMust match Vcf Header exactly``")
	var sampleName *bool = flag.Bool("sampleNames", false, "Get names of samples that appear in Vcf header")

	var samSplit *string = flag.String("sam", "", "Provide a sam alignment of the F1 hybrid to split file two sam alignments``")
	var prefix *string = flag.String("prefix", "", "Name output sam files by a defined prefix``")
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		if *sampleName && len(flag.Args()) == 1 {
			file := fileio.EasyOpen(flag.Arg(0))
			defer file.Close()
			header := vcf.ReadHeader(file)
			log.Printf("%s", vcf.PrintSampleNames(header))
		} else if *samSplit != "" {
			GoRoutinesSnpSearch(*samSplit, flag.Arg(0), parental[0], parental[1], *prefix, *f1Genome, 4)
			//SnpSearch(*samSplit, flag.Arg(0), *f1Genome, parental[0], parental[1], *prefix)
		} else {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
		}
	} else {
		input, output := flag.Arg(0), flag.Arg(1)
		switch true {
		case *filterAs:
			if len(parental) != 2 {
				log.Fatalf("Error: Must provide exactly 2 parents, found %d...\n", len(parental))
			}
			file := fileio.EasyOpen(input)
			defer file.Close()

			header := vcf.ReadHeader(file)
			sampleHash := vcf.HeaderToMaps(header)

			var parentalOne, parentalTwo, fOne int16 = sampleHash.IndexAllele[parental[0]], sampleHash.IndexAllele[parental[1]], sampleHash.IndexAllele[*f1Genome]
			writer := fileio.EasyCreate(output)
			vcf.NewWriteHeader(writer, header)
			reader := make(chan *vcf.Vcf)
			go vcf.ReadToChan(file, reader)

			defer writer.Close()
			log.SetFlags(0)
			for each := range reader {
				if vcf.ASFilter(each, parentalOne, parentalTwo, fOne) {
					//vcf.WriteVcf(writer, ReorderSampleColumns(each, int16{parentalOne, parentalTwo, fOne}))
					vcf.PrintReOrder(each, []int16{parentalOne, parentalTwo, fOne})
					vcf.WriteVcf(writer, each)
				}
			}
		}
	}
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
func SnpSearch(samfile string, genotypeVcf string, cross string, alleleOne string, alleleTwo string, prefix string) {
	var wg sync.WaitGroup

	gvcf := make(chan *vcf.Vcf)
	genotypeReader := fileio.EasyOpen(genotypeVcf)
	defer genotypeReader.Close()
	vHeader := vcf.ReadHeader(genotypeReader)
	dict := vcf.HeaderToMaps(vHeader)


	go vcf.ReadToChan(genotypeReader, gvcf)
	//children := []string{cross}
	//parents := []string{alleleOne, alleleTwo}

	//hets, homs := MapNameToIndex(dict.IndexAllele, children), MapNameToIndex(dict.IndexAllele, parents)

	snpDb := make(map[uint64]*vcf.GVcf)

	for genotype := range gvcf {
		if vcf.ASFilter(genotype, dict.IndexAllele[alleleOne], dict.IndexAllele[alleleTwo],dict.IndexAllele[cross]) {
			snpDb = vcf.GenotypeToMap(genotype, dict.FaIndex)
		}
	}
	samFile := fileio.EasyOpen(samfile)

	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	childOne := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, alleleOne))
	defer childOne.Close()
	childTwo := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, alleleTwo))
	defer childTwo.Close()

	sam.WriteHeaderToFileHandle(childOne, header)
	sam.WriteHeaderToFileHandle(childTwo, header)

	var i, parentAllele1, parentAllele2 int
	var target, query, j int64
	var ok bool
	var code uint64
	wg.Wait()

	var gV *vcf.GVcf
	for read, done := sam.NextAlignment(samFile); done != true; read, done = sam.NextAlignment(samFile) {
		parentAllele1, parentAllele2 = 0, 0
		target = read.Pos - 1
		query = 0
		for i = 0; i < len(read.Cigar); i++ {
			switch read.Cigar[i].Op {
			case 'S':
				query += read.Cigar[i].RunLength
			case 'I':
				code = vcf.ChromPosToUInt64(int(dict.FaIndex[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleTwo]) == 0 {
						parentAllele1++
					}
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleTwo]) == 0 {
						parentAllele2++
					}
				}
				query += read.Cigar[i].RunLength
			case 'D':
				code = vcf.ChromPosToUInt64(int(dict.FaIndex[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CountBase(gV.Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = vcf.ChromPosToUInt64(int(dict.FaIndex[read.RName]), int(target+j))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[dict.IndexAllele[alleleOne]].AlleleTwo]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[dict.IndexAllele[alleleTwo]].AlleleTwo]) == 0 {
							parentAllele2++
						}
					}

				}
				target += read.Cigar[i].RunLength
				query += read.Cigar[i].RunLength
			}
		}
		if parentAllele1 > parentAllele2 {
			sam.WriteAlnToFileHandle(childOne, read)
		} else if parentAllele2 > parentAllele1 {
			sam.WriteAlnToFileHandle(childTwo, read)
		} else {
			//Skip read
		}
	}
}*/

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
				code = vcf.ChromPosToUInt64(int(sampleHash.FaIndex[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[0]]].AlleleTwo]) == 0 {
						parentAllele1++
					}
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Seq[gV.Genotypes[sampleHash.IndexAllele[parents[1]]].AlleleTwo]) == 0 {
						parentAllele2++
					}
				}
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
}
