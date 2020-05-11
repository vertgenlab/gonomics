package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type Genotype struct {
	Alleles [][]dna.Base
	GVcf    []Haplotype
	//Maybe Qual?
	//I personally think it is up to the user
	//to QC input data
}

type Haplotype struct {
	One    int32
	Two    int32
	Phased bool
}

type Dictionary struct {
	Fa     map[string]int32
	HapIdx map[string]int32
}

func GenotypeToMap(vcfs []*Vcf, names map[string]int) map[uint64]*Genotype {
	mapToGVcf := make(map[uint64]*Genotype)
	var code uint64
	for _, v := range vcfs {
		code = secretCode(names[v.Chr], int(v.Pos))
		_, ok := mapToGVcf[code]
		if !ok {
			mapToGVcf[code] = vcfToGenotype(v)
		}
	}
	return mapToGVcf
}

//tmp , this functions lives in simple graph, but import cycles are not allowed...
//need to find a new package for this function
func secretCode(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

//Helper functions to convert Vcf line into a more compact version
//to accommadate a magnitude of samples.
func vcfToGenotype(v *Vcf) *Genotype {
	gVcf := &Genotype{Alleles: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), GVcf: notesToHaplotype(genotypeHelper(v.Notes))}
	return gVcf
}

func genotypeHelper(notes string) []string {
	text := strings.Split(notes, "\t")
	//var hap string
	var answer []string = make([]string, len(text))
	for i := 0; i < len(text); i++ {
		answer[i] = strings.Split(text[i], ":")[0]
	}
	return answer
}

func getAltBases(alt string) [][]dna.Base {
	words := strings.Split(alt, ",")
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

func notesToHaplotype(genotypes []string) []Haplotype {
	var answer []Haplotype = make([]Haplotype, len(genotypes))
	var alleles []string
	for i := 0; i < len(answer); i++ {
		if strings.Compare(genotypes[i], "./.") == 0 || strings.Compare(genotypes[i], ".|.") == 0 {
			answer[i] = Haplotype{One: -1, Two: -1, Phased: false}
		} else if strings.Contains(genotypes[i], "|") {
			alleles = strings.SplitN(genotypes[i], "|", 2)
			answer[i] = Haplotype{One: common.StringToInt32(alleles[0]), Two: common.StringToInt32(alleles[1]), Phased: true}
		} else {
			alleles = strings.SplitN(genotypes[i], "/", 2)
			answer[i] = Haplotype{One: common.StringToInt32(alleles[0]), Two: common.StringToInt32(alleles[1]), Phased: false}
		}
	}
	return answer
}

func HeaderToMaps(filename string) *Dictionary {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	var line, name string
	var index, hapIdx int32
	var err error
	var nextBytes []byte
	var hash *Dictionary = &Dictionary{Fa: make(map[string]int32), HapIdx: make(map[string]int32)}
	for nextBytes, err = file.Peek(1); nextBytes[0] == '#' && err == nil; nextBytes, err = file.Peek(1) {
		line, _ = fileio.EasyNextLine(file)
		if strings.HasPrefix(line, "##contig") {
			//##contig=<ID=scaffold_15,length=15191604>
			//split to get chrom/contig name
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.Fa[name]
			if !ok {
				hash.Fa[name] = index
				index++
			}
		}
		if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int32(len(words)); hapIdx++ {
				hash.HapIdx[words[hapIdx]] = hapIdx
			}
		}
		if err != nil {
			log.Fatal("There was an error reading the header line")
		}
	}
	return hash
}

//To string functions for debugging and visualizing
func genotypeToString(gVcf *Genotype) string {
	var answer string
	answer = fmt.Sprintf("0=%s,1=%s,%v", dna.BasesToString(gVcf.Alleles[0]), altBasesToString(gVcf.Alleles[1:]), haplotypesToString(gVcf))
	return answer
}

func haplotypesToString(sample *Genotype) string {
	var answer string = ""
	for i := 0; i < len(sample.GVcf); i++ {
		if sample.GVcf[i].One < 0 {
			answer += "NoData "
		} else {
			answer += fmt.Sprintf("%v=%s,%s ", sample.GVcf[i], dna.BasesToString(sample.Alleles[sample.GVcf[i].One]), dna.BasesToString(sample.Alleles[sample.GVcf[i].Two]))
		}
	}
	return answer
}

func altBasesToString(alt [][]dna.Base) string {
	var work []string = make([]string, len(alt))
	for i := 0; i < len(work); i++ {
		work[i] = dna.BasesToString(alt[i])
	}
	return strings.Join(work, ",")
}

/*
func chromPosToString(chr string, pos int64) string {
	return fmt.Sprintf("%s_%d", chr, pos)
}*/
