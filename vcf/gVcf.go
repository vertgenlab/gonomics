package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"sync"
)

//Could also call it VcfChan or HeavyVcf
type GVcf struct {
	Reader *fileio.EasyReader
	Header *VcfHeader
	Vcfs   chan *Vcf
	SyncWg *sync.WaitGroup
}

type Genotypes struct {
	Vcf
	Seq       [][]dna.Base
	Genotypes []Sample
}

type Sample struct {
	AlleleOne int16
	AlleleTwo int16
	Phased    bool
}

type SampleIdMap struct {
	FaIndex     map[string]int16
	IndexAllele map[string]int16
}

func ReadGVcf(filename string) *GVcf {
	var ans *GVcf = &GVcf{}
	var wg sync.WaitGroup
	ans.Reader = fileio.EasyOpen(filename)
	ans.Header = ReadHeader(ans.Reader)
	ans.Vcfs = make(chan *Vcf)
	ans.SyncWg = &wg
	go ReadToChan(ans.Reader, ans.Vcfs)
	return ans
}

//Uses Vcf header to create 2 hash maps 1) is the sample index that maps the which allele each sample has in Vcf 2) hash reference chromsome names to an index (used to build uint64 containing chromID and position)
func HeaderToMaps(header *VcfHeader) *SampleIdMap {
	var name string
	var index, hapIdx int16
	var hash *SampleIdMap = &SampleIdMap{FaIndex: make(map[string]int16), IndexAllele: make(map[string]int16)}
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
				hash.IndexAllele[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
}
//TODO: Solidify names
func VcfToGVcf(v *Vcf) *Genotypes {
	return &Genotypes{Vcf: *v, Seq: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), Genotypes: GetAlleleGenotype(v)}
}

func getAltBases(alt string) [][]dna.Base {
	words := strings.Split(alt, ",")
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

func GetAlleleGenotype(v *Vcf) []Sample {
	text := strings.Split(v.Notes, "\t")
	var hap string
	var alleles []string
	var answer []Sample = make([]Sample, len(text))
	for i := 0; i < len(text); i++ {
		hap = strings.Split(text[i], ":")[0]
		if strings.Compare(hap, "./.") == 0 || strings.Compare(hap, ".|.") == 0 {
			answer[i] = Sample{AlleleOne: -1, AlleleTwo: -1, Phased: false}
		} else if strings.Contains(hap, "|") {
			alleles = strings.SplitN(hap, "|", 2)
			answer[i] = Sample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true}
		} else {
			alleles = strings.SplitN(hap, "/", 2)
			answer[i] = Sample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false}
		}
	}
	return answer
}

//builds hash to Gvcfs, for a single vcf line
func GenotypeToMap(v *Vcf, names map[string]int16) map[uint64]*Genotypes {
	mapToGVcf := make(map[uint64]*Genotypes)
	return buildGenotypeMap(v, names, mapToGVcf)
}

//Parse Vcf header to quickly print sample names that appear inside Vcf
func PrintSampleNames(header *VcfHeader) string {
	var ans string = ""
	for _, line := range header.Text {
		if strings.HasPrefix(line, "#CHROM") {
			//starts at first sample column
			words := strings.Split(line, "\t")
			for i := 9; i < len(words); i++ {
				ans += fmt.Sprintf("%s\n", words[i])
			}
			return ans
		}
	}
	return ans
}

func buildGenotypeMap(v *Vcf, names map[string]int16, mapToGVcf map[uint64]*Genotypes) map[uint64]*Genotypes {
	code := ChromPosToUInt64(int(names[v.Chr]), int(v.Pos-1))
	_, ok := mapToGVcf[code]
	if !ok {
		mapToGVcf[code] = VcfToGVcf(v)
	}
	return mapToGVcf
}

func ChromPosToUInt64(chrom int, start int) uint64 {
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
	Genotypes := VcfToGVcf(v)
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, genotypeToString(Genotypes))
}

func PrintReOrder(v *Vcf, samples []int16) {
	Genotypes := VcfToGVcf(ReorderSampleColumns(v, samples))
	log.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, genotypeToString(Genotypes))
}

func vcfPrettyPrint(v *Vcf) {
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, v.Notes)
}

func genotypeToString(sample *Genotypes) string {
	var answer string = ""
	for i := 0; i < len(sample.Genotypes); i++ {
		answer += helperGenotypeToString(sample, i)
	}
	return answer
}

func helperGenotypeToString(sample *Genotypes, i int) string {
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
