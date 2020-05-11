package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

type GVcf struct {
	Alleles   [][]dna.Base
	Genotypes []Haplotype
}

type Haplotype struct {
	One    int16
	Two    int16
	Phased bool
	//GQ     uint8
}

type Dictionary struct {
	Fa     map[string]int16
	HapIdx map[string]int16
}

func GenotypeHelper(v *Vcf) []Haplotype {
	text := strings.Split(v.Notes, "\t")
	var hap string
	var alleles []string
	var answer []Haplotype = make([]Haplotype, len(text))
	for i := 0; i < len(text); i++ {
		hap = strings.Split(text[i], ":")[0]
		if strings.Compare(hap, "./.") == 0 || strings.Compare(hap, ".|.") == 0 {
			//answer[i] = Haplotype{One: -1, Two: -1, Phased: false, GQ: getGQ(v)}
			answer[i] = Haplotype{One: -1, Two: -1, Phased: false}
		} else if strings.Contains(hap, "|") {
			alleles = strings.SplitN(hap, "|", 2)
			//answer[i] = Haplotype{One: common.StringToInt16(alleles[0]), Two: common.StringToInt16(alleles[1]), Phased: true, GQ: getGQ(v)}
			answer[i] = Haplotype{One: common.StringToInt16(alleles[0]), Two: common.StringToInt16(alleles[1]), Phased: true}
		} else {
			alleles = strings.SplitN(hap, "/", 2)
			//answer[i] = Haplotype{One: common.StringToInt16(alleles[0]), Two: common.StringToInt16(alleles[1]), Phased: false, GQ: getGQ(v)}
			answer[i] = Haplotype{One: common.StringToInt16(alleles[0]), Two: common.StringToInt16(alleles[1]), Phased: false}
		}
	}
	return answer
}

func GenotypeToMap(v *Vcf, names map[string]int16, mapToGVcf map[uint64]*GVcf) map[uint64]*GVcf {
	//mapToGVcf :=
	//var code uint64
	code := secretCode(int(names[v.Chr]), int(v.Pos-1))
	_, ok := mapToGVcf[code]
	if !ok {
		mapToGVcf[code] = VcfToGenotype(v)
	}
	return mapToGVcf
}

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
func compareHaps(a Haplotype, b Haplotype) bool {
	if a.One == b.One && a.Two == b.Two {
		return true
	} else {
		return false
	}
}

func CompareAllHaps(gt []Haplotype) bool {
	for i := 0; i < len(gt)-1; i++ {
		if !compareHaps(gt[i], gt[i+1]) {
			return false
		}
	}
	return true
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
//Chr: v.Chr, Pos: int(v.Pos),
func VcfToGenotype(v *Vcf) *GVcf {
	gVcf := &GVcf{Alleles: append([][]dna.Base{dna.StringToBases(v.Ref)}, getAltBases(v.Alt)...), Genotypes: GenotypeHelper(v)}
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

func GetHeader(filename string) *VcfHeader {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	header := ReadHeader(file)
	return header
}

func HeaderToMaps(file *fileio.EasyReader) *Dictionary {
	header := ReadHeader(file)
	var name string
	var index, hapIdx int16
	var hash *Dictionary = &Dictionary{Fa: make(map[string]int16), HapIdx: make(map[string]int16)}
	for _, line := range header.Text {
		if strings.HasPrefix(line, "##contig") {
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.Fa[name]
			if !ok {
				hash.Fa[name] = index
				index++
			}
		}
		if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int16(len(words)); hapIdx++ {
				hash.HapIdx[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
}

func hapToStringHelper(sample *GVcf, i int) string {
	if sample.Genotypes[i].One < 0 {
		return "NoData "
	} else {
		if i == len(sample.Genotypes)-1 {
			return fmt.Sprintf("%d%s%d=%s%s%s", sample.Genotypes[i].One, PhasedToString(sample.Genotypes[i].Phased), sample.Genotypes[i].Two, dna.BasesToString(sample.Alleles[sample.Genotypes[i].One]), PhasedToString(sample.Genotypes[i].Phased), dna.BasesToString(sample.Alleles[sample.Genotypes[i].Two]))
		} else {
			return fmt.Sprintf("%d%s%d=%s%s%s\t", sample.Genotypes[i].One, PhasedToString(sample.Genotypes[i].Phased), sample.Genotypes[i].Two, dna.BasesToString(sample.Alleles[sample.Genotypes[i].One]), PhasedToString(sample.Genotypes[i].Phased), dna.BasesToString(sample.Alleles[sample.Genotypes[i].Two]))
			//return fmt.Sprintf("GQ=%d,%d%s%d=%s%s%s\t", sample.Genotypes[i].GQ, sample.Genotypes[i].One, PhasedToString(sample.Genotypes[i].Phased), sample.Genotypes[i].Two, dna.BasesToString(sample.Alleles[sample.Genotypes[i].One]), PhasedToString(sample.Genotypes[i].Phased), dna.BasesToString(sample.Alleles[sample.Genotypes[i].Two]))
		}
	}
}

func haplotypesToString(sample *GVcf) string {
	var answer string = ""
	for i := 0; i < len(sample.Genotypes); i++ {
		answer += hapToStringHelper(sample, i)
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

func PhasedToString(phased bool) string {
	if phased {
		return "|"
	} else {
		return "/"
	}
}
