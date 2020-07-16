package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
)

func SnpSearch(samfile string, genotypeVcf string, fOne string, parentOne string, parentTwo, prefix string) {
	reader := vcf.GoReadGVcf(genotypeVcf)
	sampleHash := vcf.HeaderToMaps(reader.Header)
	snpDb := make(map[uint64]*vcf.GVcf)
	for genotype := range reader.Vcfs {
		if vcf.ASFilter(genotype, sampleHash.GIndex[parentOne], sampleHash.GIndex[parentTwo], sampleHash.GIndex[fOne]) {
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
	var target, query, j int64
	var ok bool
	var code uint64

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
				code = vcf.ChromPosToUInt64(int(sampleHash.Fa[read.RName]), int(target))
				gV, ok = snpDb[code]
				if ok {
					if dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = vcf.ChromPosToUInt64(int(sampleHash.Fa[read.RName]), int(target+j))
					gV, ok = snpDb[code]
					if ok {
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[sampleHash.GIndex[parentOne]].AlleleTwo]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[sampleHash.GIndex[parentTwo]].AlleleTwo]) == 0 {
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
			sam.WriteAlnToFileHandle(childOne, read)
		case parentAllele2 > parentAllele1:
			sam.WriteAlnToFileHandle(childTwo, read)
		}
	}
	reader.File.Close()
}
