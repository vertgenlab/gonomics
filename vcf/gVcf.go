package vcf

import (
	"fmt"
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	//"github.com/vertgenlab/gonomics/fileio"
	"log"
	//"strconv"
	"strings"
	//"sync"
)

/*
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
	wg.Add(1)
	go ReadToChan(ans.File, ans.Vcfs, ans.SyncWg)

	go func() {
		wg.Wait()
		close(ans.Vcfs)
	}()

	return ans
}*/

type GVcf struct { //TODO: Uncommented for now, but this struct needs to be removed soon.
	Vcf
	Seq       [][]dna.Base
	Genotypes []GenomeSample
}

type SampleHash struct {
	Fa     map[string]int16
	GIndex map[string]int16
}

//TODO: Can only process short variants. Need long term solution for large structural variance.
func VcfToGvcf(v *Vcf) *GVcf {
	gVcf := &GVcf{Vcf: *v, Seq: append([][]dna.Base{dna.StringToBases(v.Ref)}, GetAltBases(v.Alt)...), Genotypes: v.Samples}
	return gVcf
}

/*
//This function has now been incorporated into ParseNotes in vcf.go
func GetAlleleGenotype(v *Vcf) []GenomeSample {
	text := strings.Split(v.Notes, "\t")
	var hap string
	var alleles []string
	var err error
	var n int64
	var answer []GenomeSample = make([]GenomeSample, len(text))
	for i := 0; i < len(text); i++ {
		hap = strings.Split(text[i], ":")[0]
		if strings.Compare(hap, "./.") == 0 || strings.Compare(hap, ".|.") == 0 {
			answer[i] = GenomeSample{AlleleOne: -1, AlleleTwo: -1, Phased: false}
		} else if strings.Contains(hap, "|") {
			alleles = strings.SplitN(hap, "|", 2)
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: true}
		} else if strings.Contains(hap, "/") {
			alleles = strings.SplitN(hap, "/", 2)
			answer[i] = GenomeSample{AlleleOne: common.StringToInt16(alleles[0]), AlleleTwo: common.StringToInt16(alleles[1]), Phased: false}
		} else {
			//Deal with single haps. There might be a better soltuion, but I think this should work.
			n, err = strconv.ParseInt(alleles[0], 10, 16)
			if err != nil && n < int64(len(text)) {
				answer[i] = GenomeSample{AlleleOne: int16(n), AlleleTwo: -1, Phased: false}
			} else {
				log.Fatalf("Error: Unexpected parsing error...\n")
			}
		}
	}
	return answer
}*/


func BuildGenotypeMap(v *Vcf, names map[string]int16, mapToVcf map[uint64]*Vcf) map[uint64]*Vcf {
	code := ChromPosToUInt64(int(names[v.Chr]), v.Pos-1)
	_, ok := mapToVcf[code]
	if !ok {
		mapToVcf[code] = v
	}
	return mapToVcf
}

/* This function is unannotated and I'm not sure what it's supposed to do. Appears to return only GQ for the first value, but TODO: this should be implemented with new VCF struct, maybe returning a slice of GQ data corresponding to each sample
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
}*/

//tmp , this functions lives in simple graph, but import cycles are not allowed...
//need to find a new package for this function
func ChromPosToUInt64(chrom int, start int) uint64 {
	var chromCode uint64 = uint64(chrom)
	chromCode = chromCode << 32
	var answer uint64 = chromCode | uint64(start)
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

func GetAltBases(words []string) [][]dna.Base {
	var answer [][]dna.Base = make([][]dna.Base, len(words))
	for i := 0; i < len(words); i++ {
		answer[i] = dna.StringToBases(words[i])
	}
	return answer
}

func AltBasesToStrings(alt [][]dna.Base) []string {
	var work []string = make([]string, len(alt))
	for i := 0; i < len(work); i++ {
		work[i] = dna.BasesToString(alt[i])
	}
	return work
}

func PhasedToString(phased bool) string {
	if phased {
		return "|"
	} else {
		return "/"
	}
}

//ReorderSampleColumns reorganizes the Samples slice based on a samples []int16 specification list.
func ReorderSampleColumns(input *Vcf, samples []int16) *Vcf {
	outSamples := make([]GenomeSample, 0, len(samples))
	for i := 0; i < len(samples); i++ {
		outSamples = append(outSamples, input.Samples[samples[i]])
	}
	input.Samples = outSamples
	return input
}

func PrintReOrder(v *Vcf, samples []int16) {
	vReorder := ReorderSampleColumns(v, samples)
	log.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, SamplesToString(vReorder.Samples))
}

func SamplesToString(sample []GenomeSample) string {
	var answer string = ""
	for i := 0; i < len(sample); i++ {
		answer += HelperSamplesToString(sample, i)
	}
	return answer
}

//helperGenotypeToStringNew uses just an array of GenomeSample structs to write to a string for simple gVCFs with just the allele info in notes.
func HelperSamplesToString(sample []GenomeSample, i int) string {
	var answer string
	if sample[i].AlleleOne < 0 {
		answer = "noData"
	} else {
		if i == len(sample)-1 {
			answer = fmt.Sprintf("%d%s%d", sample[i].AlleleOne, PhasedToString(sample[i].Phased), sample[i].AlleleTwo)
		} else {
			answer = fmt.Sprintf("%d%s%d\t", sample[i].AlleleOne, PhasedToString(sample[i].Phased), sample[i].AlleleTwo)
		}
	}
	if len(sample[i].FormatData) > 0 {
		answer = answer + strings.Join(sample[i].FormatData, ":")
	}
	return answer
}

func vcfPrettyPrint(v *Vcf) {
	fmt.Printf("%s\t%d\t%s\t%s\t%s\n", v.Chr, v.Pos, v.Ref, v.Alt, SamplesToString(v.Samples))
}
