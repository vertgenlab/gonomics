package vcf

import (
	"fmt"
	//"strconv".
	"strings"

	//"github.com/vertgenlab/gonomics/common".
	"github.com/vertgenlab/gonomics/dna"
)

type GVcf struct { //TODO: Uncommented for now, but this struct needs to be removed soon.
	Vcf
	Seq       [][]dna.Base
	Genotypes []Sample
}

type SampleHash struct {
	Fa     map[string]int16
	GIndex map[string]int16
}

func BuildGenotypeMap(v Vcf, names map[string]int16, mapToVcf map[uint64]Vcf) map[uint64]Vcf {
	code := ChromPosToUInt64(int(names[v.Chr]), v.Pos-1)
	_, ok := mapToVcf[code]
	if !ok {
		mapToVcf[code] = v
	}
	return mapToVcf
}

// tmp , this functions lives in simple graph, but import cycles are not allowed...
// need to find a new package for this function.
func ChromPosToUInt64(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
	return answer
}

// Parse Vcf header to quickly print sample names that appear inside Vcf.
func PrintSampleNames(header Header) string {
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

func GetAltBases(words []string) [][]dna.Base {
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

func PhasedToString(phased bool) string {
	if phased {
		return "|"
	} else {
		return "/"
	}
}

// ReorderSampleColumns reorganizes the Samples slice based on a samples []int16 specification list.
func ReorderSampleColumns(input Vcf, samples []int16) Vcf {
	outSamples := make([]Sample, 0, len(samples))
	for i := 0; i < len(samples); i++ {
		outSamples = append(outSamples, input.Samples[samples[i]])
	}
	input.Samples = outSamples
	return input
}

func SamplesToString(sample []Sample) string {
	var answer string = ""
	for i := 0; i < len(sample); i++ {
		if i > 0 {
			answer += "\t" + sampleToString(sample[i])
		} else {
			answer += sampleToString(sample[i])
		}
	}
	return answer
}

// sampleToString uses just an array of Sample structs to write to a string for simple gVCFs with just the allele info in notes.
func sampleToString(s Sample) string {
	var answer string
	if s.FormatData == nil {
		return "."
	}
	if s.Alleles == nil {
		answer = "."
	} else {
		answer += fmt.Sprintf("%d", s.Alleles[0])
		for i := 1; i < len(s.Phase); i++ {
			answer += fmt.Sprintf("%s%d", PhasedToString(s.Phase[i]), s.Alleles[i])
		}
	}
	if len(s.FormatData) > 0 {
		answer = answer + strings.Join(s.FormatData, ":")
	}

	return answer
}
