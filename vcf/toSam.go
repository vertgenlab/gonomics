package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	//"github.com/vertgenlab/gonomics/common"
	"fmt"
	//"log"
	"os"
	"sync"
	"strings"
)

func SnpSearch(samfile string, genotypeVcf string, cross string, alleleOne string, alleleTwo, prefix string) {
	dict := HeaderToMaps(genotypeVcf)

	var wg sync.WaitGroup
	gvcf := make(chan *Vcf)
	go ReadToChan(genotypeVcf, gvcf)
	Aa:= strings.Split(cross, ",")
	aa := []string{alleleOne, alleleTwo}
	// ReadFilterList(sampleSheet, dict.HapIdx)
	index1, index2 := MapNameToIndex(dict.HapIdx, Aa), MapNameToIndex(dict.HapIdx, aa)

	samFile := fileio.EasyOpen(samfile)
	defer samFile.Close()
	header := sam.ReadHeader(samFile)

	AA, _ := os.Create(fmt.Sprintf("%s.%s.SNPs.sam", prefix,aa[0]))
	defer AA.Close()
	bb, _ := os.Create(fmt.Sprintf("%s.%s.SNPs.sam", prefix,aa[1]))
	defer bb.Close()

	sam.WriteHeaderToFileHandle(AA, header)
	sam.WriteHeaderToFileHandle(bb, header)
	wg.Add(1)

	snpDb := GenotypeToMap(gvcf, dict.Fa, index1, index2, &wg)
	var j, refCount, altCount int
	var target, query, k int64
	var ok bool
	var code uint64

	wg.Wait()
	for read, done := sam.NextAlignment(samFile); done != true; read, done = sam.NextAlignment(samFile) {
		refCount, altCount = 0, 0
		target = read.Pos - 1
		query = 0
		for j = 0; j < len(read.Cigar); j++ {
			switch read.Cigar[j].Op {
			case 'S':
				query += read.Cigar[j].RunLength
			case 'I':
				query += read.Cigar[j].RunLength
			case 'D':
				code = secretCode(int(dict.Fa[read.RName]), int(target))
				_, ok = snpDb[code]
				if ok {
					if dna.Gap == snpDb[code].Alleles[snpDb[code].Genotypes[dict.HapIdx[aa[0]]].One][0] {
						refCount++
					}
				}
				target += read.Cigar[j].RunLength
			case 'M':
				for k = 0; k < read.Cigar[j].RunLength; k++ {
					code = secretCode(int(dict.Fa[read.RName]), int(target+k))
					_, ok = snpDb[code]
					if ok {
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+k]}, snpDb[code].Alleles[snpDb[code].Genotypes[dict.HapIdx[aa[0]]].One]) == 0 {
							refCount++
						}
						if dna.CompareSeqsIgnoreCase([]dna.Base{read.Seq[query+k]}, snpDb[code].Alleles[snpDb[code].Genotypes[dict.HapIdx[aa[1]]].One]) == 0 {
							altCount++
						}
					}
				}
				target += read.Cigar[j].RunLength
				query += read.Cigar[j].RunLength
			}
		}
		if refCount > altCount {
			sam.WriteAlnToFileHandle(AA, read)
		} else if altCount > refCount {
			sam.WriteAlnToFileHandle(bb, read)
		} else {
			//Skip read
		}
	}
}
