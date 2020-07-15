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

type Reader struct {
	File   *fileio.EasyReader
	Header *VcfHeader
	Vcfs   chan *Vcf
	SyncWg *sync.WaitGroup
}

func GoReadGVcf(filename string) *Reader {
	var ans *Reader = &Reader{}
	var wg sync.WaitGroup
	ans.File = fileio.EasyOpen(filename)
	ans.Header = ReadHeader(ans.File)
	ans.Vcfs = make(chan *Vcf)
	ans.SyncWg = &wg
	go ReadToChan(ans.File, ans.Vcfs)
	return ans
}

type GVcf struct {
	Vcf
	Seq       [][]dna.Base
	Genotypes []GenomeSample
}

type GenomeSample struct {
	AlleleOne int16
	AlleleTwo int16
	Phased    bool
}

type SampleHash struct {
	Fa     map[string]int16
	GIndex map[string]int16
}

func GetAlleleGenotype(v *Vcf) []GenomeSample {
	text := strings.Split(v.Notes, "\t")
	var hap string
	var alleles []string
	var answer []GenomeSample = make([]GenomeSample, len(text))
	for i := 0; i < len(text); i++ {
		hap = strings.Split(text[i], ":")[0]
		if strings.Compare(hap, "./.") == 0 || strings.Compare(hap, ".|.") == 0 {
			//answer[i] = GenomeSample{AlleleOne: -1, AlleleTwo: -1, Phased: false, GQ: getGQ(v)}
			answer[i] = GenomeSample{AlleleOne: -1, AlleleTwo: -1, Phased: false}
		} else if strings.Contains(hap, "|") {
			alleles = strings.SplitN(hap, "|", 2)
			//answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true, GQ: getGQ(v)}
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true}
		} else {
			alleles = strings.SplitN(hap, "/", 2)
			//answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false, GQ: getGQ(v)}
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false}
		}
	}
	return answer
}

func BuildGenotypeMap(v *Vcf, names map[string]int16, mapToGVcf map[uint64]*GVcf) map[uint64]*GVcf {
	code := ChromPosToUInt64(int(names[v.Chr]), int(v.Pos-1))
	_, ok := mapToGVcf[code]
	if !ok {
		mapToGVcf[code] = VcfToGvcf(v)
	}
	return mapToGVcf
}

/*
//builds hash to Gvcfs, for a single vcf line
func GenotypeToMap(v *Vcf, names map[string]int16) map[uint64]*GVcf {
	mapToGVcf := make(map[uint64]*GVcf)
	return buildGenotypeMap(v, names, mapToGVcf)
}

func buildGenotypeMap(v *Vcf, names map[string]int16, mapToGVcf map[uint64]*GVcf) map[uint64]*GVcf {
	code := ChromPosToUInt64(int(names[v.Chr]), int(v.Pos-1))
	_, ok := mapToGVcf[code]
	if !ok {
		mapToGVcf[code] = VcfToGVcf(v)
	}
	return mapToGVcf
}*/

func getGQ(v *Vcf) uint8 {
	var answer uint8 = 0
	if strings.Contains(v.Format, "GQ") {
		stats := strings.Split(v.Format, ":")
		for i := 0; i < len(stats); i++ {
			if strings.Compare(stats[i], "GQ") == 0 {
				value := strings.Split(v.Notes, ":")
				if strings.Contains(value[i], ".") || strings.Contains(value[i], ",") {
					answer = 0
				} else {

					answer = common.StringToUint8(value[i])
				}

			}

		}
	}
	return answer
}
func compareHaps(a GenomeSample, b GenomeSample) bool {
	if a.AlleleOne == b.AlleleOne && a.AlleleTwo == b.AlleleTwo {
		return true
	} else {
		return false
	}
}

func CompareAllHaps(gt []GenomeSample) bool {
	for i := 0; i < len(gt)-1; i++ {
		if !compareHaps(gt[i], gt[i+1]) {
			return false
		}
	}
	return true
}

//tmp , this functions lives in simple graph, but import cycles are not allowed...
//need to find a new package for this function
func ChromPosToUInt64(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

//Helper functions to convert Vcf line into a more compact version
//to accommadate a magnitude of samples.
//Chr: v.Chr, Pos: int(v.Pos),
func VcfToGvcf(v *Vcf) *GVcf {
	gVcf := &GVcf{Vcf: *v, Seq: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), Genotypes: GetAlleleGenotype(v)}
	return gVcf
}

//func VcfToGVcf(v *Vcf) *Genotypes {
//	return &Genotypes{Vcf: *v, Seq: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), Genotypes: GetAlleleGenotype(v)}
//}

func getAltBases(alt string) [][]dna.Base {
	words := strings.Split(alt, ",")
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

//Uses Vcf header to create 2 hash maps 1) is the sample index that maps the which allele each sample has in Vcf 2) hash reference chromsome names to an index (used to build uint64 containing chromID and position)
func HeaderToMaps(header *VcfHeader) *SampleHash {
	var name string
	var index, hapIdx int16
	var hash *SampleHash = &SampleHash{Fa: make(map[string]int16), GIndex: make(map[string]int16)}
	for _, line := range header.Text {
		if strings.HasPrefix(line, "##contig") {
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.Fa[name]
			if !ok {
				hash.Fa[name] = index
				index++
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int16(len(words)); hapIdx++ {
				hash.GIndex[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
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

func PrintReOrder(v *Vcf, samples []int16) {
	Genotypes := VcfToGvcf(ReorderSampleColumns(v, samples))
	log.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, genotypeToString(Genotypes))
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

func vcfPrettyPrint(v *Vcf) {
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, v.Notes)
}
