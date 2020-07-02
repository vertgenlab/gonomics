package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"strings"
)

type GVcf struct {
	Seq       [][]dna.Base
	Genotypes []Genotype
}

type Genotype struct {
	AlleleOne int16
	AlleleTwo int16
	Phased    bool
}

type SampleIdMap struct {
	FaIndex     map[string]int16
	IndexAlelle map[string]int16
}

func HeaderToMaps(header *VcfHeader) *SampleIdMap {
	var name string
	var index, hapIdx int16
	var hash *SampleIdMap = &SampleIdMap{FaIndex: make(map[string]int16), IndexAlelle: make(map[string]int16)}
	for _, line := range header.Text {
		if strings.HasPrefix(line, "##contig") {
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.FaIndex[name]
			if !ok {
				hash.FaIndex[name] = index
				index++
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int16(len(words)); hapIdx++ {
				hash.IndexAlelle[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
}

func VcfToGenotype(v *Vcf) *GVcf {
	gVcf := &GVcf{Seq: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), Genotypes: GetAlleleGenotype(v)}
	return gVcf
}

func getAltBases(alt string) [][]dna.Base {
	words := strings.Split(alt, ",")
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

func GetAlleleGenotype(v *Vcf) []Genotype {
	text := strings.Split(v.Notes, "\t")
	var hap string
	var alleles []string
	var answer []Genotype = make([]Genotype, len(text))
	for i := 0; i < len(text); i++ {
		hap = strings.Split(text[i], ":")[0]
		if strings.Compare(hap, "./.") == 0 || strings.Compare(hap, ".|.") == 0 {
			answer[i] = Genotype{AlleleOne: -1, AlleleTwo: -1, Phased: false}
		} else if strings.Contains(hap, "|") {
			alleles = strings.SplitN(hap, "|", 2)
			answer[i] = Genotype{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true}
		} else {
			alleles = strings.SplitN(hap, "/", 2)
			answer[i] = Genotype{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false}
		}
	}
	return answer
}

func GenotypeToMap(v *Vcf, names map[string]int16) map[uint64]*GVcf {
	mapToGVcf := make(map[uint64]*GVcf)
	return buildGenotypeMap(v, names, mapToGVcf)
}

func buildGenotypeMap(v *Vcf, names map[string]int16, mapToGVcf map[uint64]*GVcf) map[uint64]*GVcf {
	code := chromPosToUInt64(int(names[v.Chr]), int(v.Pos-1))
	_, ok := mapToGVcf[code]
	if !ok {
		mapToGVcf[code] = VcfToGenotype(v)
	}
	return mapToGVcf
}

func chromPosToUInt64(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

func ReorderSampleColumns(input *Vcf, samples []int16) *Vcf {
	columnData := make([]string, 0, len(samples))
	words := strings.Split(input.Notes, "\t")
	var i int
	for i = 0; i < len(samples); i++ {
		columnData = append(columnData, words[samples[i]])
	}
	input.Notes = strings.Join(columnData, "\t")
	return input
}

func ViewGenotypeVcf(v *Vcf) {
	gVcf := VcfToGenotype(v)
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, genotypeToString(gVcf))
}

func PrintReOrder(v *Vcf, samples []int16) {
	gVcf := VcfToGenotype(ReorderSampleColumns(v, samples))
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, genotypeToString(gVcf))
}

func vcfPrettyPrint(v *Vcf) {
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, v.Notes)
}

func genotypeToString(sample *GVcf) string {
	var answer string = ""
	for i := 0; i < len(sample.Genotypes); i++ {
		answer += helperGenotypeToString(sample, i)
	}
	return answer
}

func helperGenotypeToString(sample *GVcf, i int) string {
	if sample.Genotypes[i].AlleleOne < 0 {
		return "NoData "
	} else {
		if i == len(sample.Genotypes)-1 {
			return fmt.Sprintf("%d%s%d=%s%s%s", sample.Genotypes[i].AlleleOne, PhasedToString(sample.Genotypes[i].Phased), sample.Genotypes[i].AlleleTwo, dna.BasesToString(sample.Seq[sample.Genotypes[i].AlleleOne]), PhasedToString(sample.Genotypes[i].Phased), dna.BasesToString(sample.Seq[sample.Genotypes[i].AlleleTwo]))
		} else {
			return fmt.Sprintf("%d%s%d=%s%s%s\t", sample.Genotypes[i].AlleleOne, PhasedToString(sample.Genotypes[i].Phased), sample.Genotypes[i].AlleleTwo, dna.BasesToString(sample.Seq[sample.Genotypes[i].AlleleOne]), PhasedToString(sample.Genotypes[i].Phased), dna.BasesToString(sample.Seq[sample.Genotypes[i].AlleleTwo]))
		}
	}
}

func altBasesToString(alt [][]dna.Base) string {
	var work []string = make([]string, len(alt))
	for i := 0; i < len(work); i++ {
		work[i] = dna.BasesToString(alt[i])
	}
	return strings.Join(work, ",")
}

func PhasedToString(phased bool) string {
	if phased {
		return "|"
	} else {
		return "/"
	}
}
