package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	//"strings"
	//"sync"
)

func SnpSearch(samfile string, genotypeVcf string, fOne string, parentOne string, parentTwo, prefix string) {
	reader := ReadGVcf(genotypeVcf)
	dict := HeaderToMaps(reader.Header)
	snpDb := make(map[uint64]*GVcf)
	for genotype := range reader.Vcfs {
		if ASFilter(genotype, dict.GIndex[parentOne], dict.GIndex[parentTwo], dict.GIndex[fOne]) {
			GenotypeToMap(genotype, dict.Fa, snpDb)
		}
	}
	samFile := fileio.EasyOpen(samfile)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	childOne := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentOne))
	defer childOne.Close()
	childTwo := fileio.EasyCreate(fmt.Sprintf("%s.%s.SNPs.sam", prefix, parentTwo))
	defer childTwo.Close()

	sam.WriteHeaderToFileHandle(childOne, header)
	sam.WriteHeaderToFileHandle(childTwo, header)

	var i, parentAllele1, parentAllele2 int
	var target, query, j int64
	var ok bool
	var code uint64
	//wg.Wait()

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
				//code = ChromPosToUInt64(int(dict.Fa[read.RName]), int(target))
				//_, ok = snpDb[code]
				//if ok {
				//	gV = snpDb[code]
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.GIndex[parentOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.GIndex[parentOne]].AlleleTwo]) == 0 {
				//		parentAllele1++
				//	}
				//	if dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.GIndex[parentTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase(read.Seq[query:query+read.Cigar[i].RunLength], gV.Alleles[gV.Genotypes[dict.GIndex[parentTwo]].AlleleTwo]) == 0 {
				//		parentAllele2++
				//	}
				//}
				query += read.Cigar[i].RunLength
			case 'D':
				code = ChromPosToUInt64(int(dict.Fa[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					gV = snpDb[code]
					if dna.CountBase(gV.Seq[gV.Genotypes[dict.GIndex[parentOne]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[dict.GIndex[parentOne]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
					if dna.CountBase(gV.Seq[gV.Genotypes[dict.GIndex[parentTwo]].AlleleOne], dna.Gap) == int(read.Cigar[i].RunLength) && dna.CountBase(gV.Seq[gV.Genotypes[dict.GIndex[parentTwo]].AlleleTwo], dna.Gap) == int(read.Cigar[i].RunLength) {
						parentAllele1++
					}
				}
				target += read.Cigar[i].RunLength
			case 'M':
				for j = 0; j < read.Cigar[i].RunLength; j++ {
					code = ChromPosToUInt64(int(dict.Fa[read.RName]), int(target+j))
					_, ok = snpDb[code]
					if ok {
						gV = snpDb[code]
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[dict.GIndex[parentOne]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[dict.GIndex[parentOne]].AlleleTwo]) == 0 {
							parentAllele1++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, gV.Seq[gV.Genotypes[dict.GIndex[parentTwo]].AlleleOne]) == 0 && dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+j]}, snpDb[code].Seq[gV.Genotypes[dict.GIndex[parentTwo]].AlleleTwo]) == 0 {
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
	reader.File.Close()
}
