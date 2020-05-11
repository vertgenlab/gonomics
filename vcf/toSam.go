package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
	"sync"
)

func SnpSearch(samfile string, genotypeVcf string, cross string, alleleOne string, alleleTwo, prefix string) {

	var wg sync.WaitGroup

	gvcf := make(chan *Vcf)
	genotypeReader := fileio.EasyOpen(genotypeVcf)
	defer genotypeReader.Close()
	dict := HeaderToMaps(genotypeReader)

	go ReadToChan(genotypeReader, gvcf)
	children := strings.Split(cross, ",")
	parents := []string{alleleOne, alleleTwo}

	hets, homs := MapNameToIndex(dict.HapIdx, children), MapNameToIndex(dict.HapIdx, parents)

	snpDb := make(map[uint64]*GVcf)

	for genotype := range gvcf {
		if ASFilter(GenotypeHelper(genotype), hets, homs) {
			GenotypeToMap(genotype, dict.Fa, snpDb)
		}
	}
	samFile := fileio.EasyOpen(samfile)

	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	childOne := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parents[0]))
	defer childOne.Close()
	childTwo := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parents[1]))
	defer childTwo.Close()

	sam.WriteHeaderToFileHandle(childOne, header)
	sam.WriteHeaderToFileHandle(childTwo, header)

	var i, parentAllele1, parentAllele2 int
	var target, query, j int64
	var ok bool
	var code uint64
	wg.Wait()

	var gV *GVcf
	for read, done := sam.NextAlignment(samFile); done != true; read, done = sam.NextAlignment(samFile) {
		parentAllele1, parentAllele2 = 0, 0
		target = read.Pos - 1
		query = 0
		for i = 0; i < len(read.Cigar); i++ {
			switch read.Cigar[i].Op {
			case 'S':
				query += read.Cigar[i].RunLength
			case 'I':
				code = secretCode(int(dict.Fa[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two]) == 0 {
						parentAllele1++
					}
					if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two]) == 0 {
						parentAllele2++
					}
				}
				query += read.Cigar[i].RunLength
			case 'D':
				code = secretCode(int(dict.Fa[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = secretCode(int(dict.Fa[read.RName]), int(target+j))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two]) == 0 {
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
}

func GoRoutinesSnpSearch(samfile string, genotypeVcf string, cross string, alleleOne string, alleleTwo, prefix string) {

	var wg sync.WaitGroup
	gvcf := make(chan *Vcf)

	genotypeReader := fileio.EasyOpen(genotypeVcf)
	defer genotypeReader.Close()
	dict := HeaderToMaps(genotypeReader)

	go ReadToChan(genotypeReader, gvcf)

	children := strings.Split(cross, ",")
	parents := []string{alleleOne, alleleTwo}
	// ReadFilterList(sampleSheet, dict.HapIdx)
	hets, homs := MapNameToIndex(dict.HapIdx, children), MapNameToIndex(dict.HapIdx, parents)

	snpDb := make(map[uint64]*GVcf)

	for genotype := range gvcf {
		if ASFilter(GenotypeHelper(genotype), hets, homs) {
			GenotypeToMap(genotype, dict.Fa, snpDb)
		}
	}
	samFile := fileio.EasyOpen(samfile)

	defer samFile.Close()
	header := sam.ReadHeader(samFile)
	//childOne := fileio.EasyCreate()
	//defer childOne.Close()
	//childTwo := fileio.EasyCreate()
	//defer childTwo.Close()

	//sam.WriteHeaderToFileHandle(childOne, header)
	//sam.WriteHeaderToFileHandle(childTwo, header)

	samReader := make(chan *sam.SamAln)
	go sam.ReadToChan(samFile, samReader)
	wg.Wait()
	var wgReader, wgWriter sync.WaitGroup
	childOne := make(chan *sam.SamAln)
	childTwo := make(chan *sam.SamAln)
	for i := 0; i < 8; i++ {
		wgReader.Add(1)
		go snpAnalysis(snpDb, dict, parents, samReader, childOne, childTwo, &wgReader)
	}

	wgWriter.Add(2)
	go sam.SamChanToFile(childOne, fmt.Sprintf("%s.%s.SNPs.sam", prefix, parents[0]), header, &wgWriter)
	go sam.SamChanToFile(childTwo, fmt.Sprintf("%s.%s.SNPs.sam", prefix, parents[1]), header, &wgWriter)
	wgReader.Wait()
	close(childOne)
	close(childOne)
	wgWriter.Wait()
}

func snpAnalysis(snpDb map[uint64]*GVcf, dict *Dictionary, parents []string, samReader <-chan *sam.SamAln, childOne chan<- *sam.SamAln, childTwo chan<- *sam.SamAln, wg *sync.WaitGroup) {
	//for read, done := sam.NextAlignment(samFile); done != true; read, done = sam.NextAlignment(samFile) {

	for read := range samReader {
		parentAllele1, parentAllele2 := 0, 0
		var target int64 = read.Pos - 1
		var query int64 = 0
		var code uint64
		var ok bool
		var gV *GVcf
		for i := 0; i < len(read.Cigar); i++ {
			switch read.Cigar[i].Op {
			case 'S':
				query += read.Cigar[i].RunLength
			case 'I':
				/*	code = secretCode(int(dict.Fa[read.RName]), int(target))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two]) == 0 {
							parentAllele2++
						}
					}*/
				query += read.Cigar[i].RunLength
			case 'D':
				code = secretCode(int(dict.Fa[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				var j int64
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = secretCode(int(dict.Fa[read.RName]), int(target+j))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].One]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Alleles[gV.Genotypes[dict.HapIdx[parents[0]]].Two]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].One]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Alleles[gV.Genotypes[dict.HapIdx[parents[1]]].Two]) == 0 {
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
