package main

import (
	"fmt"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
)

func ASFilter(v vcf.Vcf, parentOne int16, parentTwo int16, F1 int16) bool {
	if vcf.IsHomozygous(v.Samples[parentOne]) && vcf.IsHomozygous(v.Samples[parentTwo]) && vcf.IsHeterozygous(v.Samples[F1]) && v.Samples[parentOne].Alleles[0] != v.Samples[parentTwo].Alleles[1] {
		return true
	} else {
		return false
	}
}

func SnpSearch(samfile string, genotypeVcf string, fOne string, parentOne string, parentTwo, prefix string) {
	reader, vcfHeader := vcf.GoReadToChan(genotypeVcf)
	sampleHash := vcf.HeaderToMaps(vcfHeader)
	snpDb := make(map[uint64]vcf.Vcf)
	for genotype := range reader {
		if ASFilter(genotype, sampleHash.GIndex[parentOne], sampleHash.GIndex[parentTwo], sampleHash.GIndex[fOne]) {
			vcf.BuildGenotypeMap(genotype, sampleHash.Fa, snpDb)
		}
	}
	samFile := fileio.EasyOpen(samfile)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	childOne, childTwo := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentOne)), fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentTwo))
	defer childOne.Close()
	defer childTwo.Close()

	sam.WriteHeaderToFileHandle(childOne, header)
	sam.WriteHeaderToFileHandle(childTwo, header)

	var i, parentAllele1, parentAllele2 int
	var target, query, j int
	var ok bool
	var code uint64

	var gV vcf.Vcf
	var alleles [][]dna.Base
	for read, done := sam.ReadNext(samFile); done != true; read, done = sam.ReadNext(samFile) {
		parentAllele1, parentAllele2 = 0, 0
		target = int(read.Pos - 1)
		query = 0
		alleles = vcf.GetAltBases(strings.Split(fmt.Sprintf("%s,%s", gV.Ref, gV.Alt), ","))
		for i = 0; i < len(read.Cigar); i++ {
			switch read.Cigar[i].Op {
			case 'S':
				query += read.Cigar[i].RunLength
			case 'I':
				//TODO: Figure out how to take insertions into account. This algorithm below should work in theory, but there is a case I can't figure out
				//code = ChromPosToUInt64(int(sampleHash.Fa[read.RName]), int(target))
				//_, ok = snpDb[code]
				//if ok {
				//	gV = snpDb[code]
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleTwo]) == 0 {
				//		parentAllele1++
				//	}
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleTwo]) == 0 {
				//		parentAllele2++
				//	}
				//}
				query += read.Cigar[i].RunLength
			case 'D':
				code = vcf.ChromPosToUInt64(int(sampleHash.Fa[read.RName]), target)
				gV, ok = snpDb[code]
				if ok {
					if dna.CountBase(alleles[gV.Samples[sampleHash.GIndex[parentOne]].Alleles[0]], dna.Gap) == read.Cigar[i].RunLength && dna.CountBase(alleles[gV.Samples[sampleHash.GIndex[parentOne]].Alleles[1]], dna.Gap) == read.Cigar[i].RunLength {
						parentAllele1++
					}
					if dna.CountBase(alleles[gV.Samples[sampleHash.GIndex[parentTwo]].Alleles[0]], dna.Gap) == read.Cigar[i].RunLength && dna.CountBase(alleles[gV.Samples[sampleHash.GIndex[parentTwo]].Alleles[1]], dna.Gap) == read.Cigar[i].RunLength {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = vcf.ChromPosToUInt64(int(sampleHash.Fa[read.RName]), target+j)
					gV, ok = snpDb[code]
					if ok {
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, alleles[gV.Samples[sampleHash.GIndex[parentOne]].Alleles[0]]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, alleles[snpDb[code].Samples[sampleHash.GIndex[parentOne]].Alleles[1]]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, alleles[gV.Samples[sampleHash.GIndex[parentTwo]].Alleles[0]]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, alleles[snpDb[code].Samples[sampleHash.GIndex[parentTwo]].Alleles[1]]) == 0 {
							parentAllele2++
						}
					}
				}
				target += read.Cigar[i].RunLength
				query += read.Cigar[i].RunLength
			}
		}
		switch true {
		case parentAllele1 > parentAllele2:
			sam.WriteToFileHandle(childOne, read)
		case parentAllele2 > parentAllele1:
			sam.WriteToFileHandle(childTwo, read)
		}
	}
}
