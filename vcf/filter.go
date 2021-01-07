package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	//"log"
	"strings"
)

//Filter returns true if a Vcf passes a set of filter criteria, false otherwise.
func Filter(v *Vcf, chrom string, minPos int, maxPos int, ref string, alt []string, minQual float64) bool {
	if !FilterRange(v, minPos, maxPos) {
		return false
	}
	if !FilterChrom(v, chrom) {
		return false
	}
	if !FilterRef(v, ref) {
		return false
	}
	if !FilterAlt(v, alt) {
		return false
	}
	if !FilterQual(v, minQual) {
		return false
	}
	return true
}

//FilterQual returns true if a Vcf quality is above an input minQual score, false otherwise.
func FilterQual(v *Vcf, minQual float64) bool {
	if v.Qual < minQual {
		return false
	}
	return true
}

//FilterAlt returns true if the Alt field of a Vcf entry matches a desired input Alt field, false otherwise. Order sensitive.
func FilterAlt(v *Vcf, alt []string) bool {
	if len(alt) > 0 && CompareAlt(v.Alt, alt) != 0 {
		return false
	}
	return true
}

//FilterRef returns true if the Ref field of a Vcf record matches an input string, false otherwise.
func FilterRef(v *Vcf, ref string) bool {
	if len(ref) > 0 && strings.Compare(v.Ref, ref) != 0 {
		return false
	}
	return true
}

//FilterRange returns true if a Vcf position lies between an input minimum and maximum position (inclusive for both min and max), false otherwise.
func FilterRange(v *Vcf, minPos int, maxPos int) bool {
	if v.Pos < minPos || v.Pos > maxPos {
		return false
	}
	return true
}

//FilterChrom returns true if the Chrom field of a Vcf record matches an input string, false otherwise.
func FilterChrom(v *Vcf, chrom string) bool {
	if chrom != "" && v.Chr != chrom {
		return false
	}
	return true
}

//TODO: This is re-implemented andf optimized on line 169. Once I can confirm the functions behave the same way, this will be removed.
func FilterAxtVcf(vcfs []*Vcf, fa []*fasta.Fasta) []*Vcf {
	split := VcfSplit(vcfs, fa)
	var answer []*Vcf
	var i, j int
	var ref []dna.Base
	var alt [][]dna.Base
	var noN bool
	for i = 0; i < len(split); i++ {
		encountered := make(map[int]bool)
		for j = 0; j < len(split[i]); j++ {
			if encountered[split[i][j].Pos] == true {
				//do not add
			} else {
				encountered[split[i][j].Pos] = true
				ref = dna.StringToBases(split[i][j].Ref)
				alt = GetAltBases(split[i][j].Alt)
				noN = true
				if dna.CountBaseInterval(ref, dna.N, 0, len(ref)) == 0 {
					for k := 0; k < len(alt); k++ {
						if dna.CountBaseInterval(alt[k], dna.N, 0, len(alt)) != 0 {
							noN = false
						}
					}
					if noN {
						answer = append(answer, split[i][j])
					}
				}
			}
		}
	}
	Sort(answer)
	return answer
}

//FilterNs removes all records from a slice of Vcfs that contain Ns.
func FilterNs(vcfs []*Vcf) []*Vcf {
	var answer []*Vcf
	var noN bool
	for i := 0; i < len(vcfs); i++ {
		noN = true
		if dna.CountBase(dna.StringToBases(vcfs[i].Ref), dna.N) > 0 {
			noN = false
		}
		for j := 0; j < len(vcfs[i].Alt); j++ {
			if dna.CountBase(dna.StringToBases(vcfs[i].Alt[j]), dna.N) > 0 {
				noN = false
			}
		}
		if noN {
			answer = append(answer, vcfs[i])
		}
	}
	return answer
}

func ASFilter(v *Vcf, parentOne int16, parentTwo int16, F1 int16) bool {
	if IsHomozygous(v.Samples[parentOne]) && IsHomozygous(v.Samples[parentTwo]) && IsHeterozygous(v.Samples[F1]) && v.Samples[parentOne].AlleleOne != v.Samples[parentTwo].AlleleOne {
		return true
	} else {
		return false
	}
}

func mergeSimilarVcf(a *Vcf, b *Vcf) *Vcf {
	mergeRecord := &Vcf{Chr: a.Chr, Pos: a.Pos, Id: a.Id, Ref: "", Qual: a.Qual, Filter: "Merged:SNP:INDEL", Info: a.Info, Format: a.Format, Samples: a.Samples}
	if len(a.Ref) < len(b.Ref) {
		mergeRecord.Ref += b.Ref
	} else {
		mergeRecord.Ref += a.Ref
	}
	if len(a.Alt) < len(b.Alt) {
		mergeRecord.Alt = b.Alt
	} else {
		mergeRecord.Alt = a.Alt
	}
	return mergeRecord
}

func LowQual(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return true
	} else {
		return false
	}
}

func IsHeterozygous(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return false
	}
	if genome.AlleleOne != genome.AlleleTwo {
		return true
	}
	if genome.AlleleOne == genome.AlleleTwo {
		return false
	}
	return false
}

func IsHomozygous(genome GenomeSample) bool {
	if genome.AlleleOne < 0 || genome.AlleleTwo < 0 {
		return false
	}
	if genome.AlleleOne == genome.AlleleTwo {
		return true
	}
	if genome.AlleleOne != genome.AlleleTwo {
		return false
	}
	return false
}

func getListIndex(header *VcfHeader, list []string) []int16 {
	sampleHash := HeaderToMaps(header)
	var listIndex []int16 = make([]int16, len(list))
	for i := 0; i < len(listIndex); i++ {
		//look up alt allele index belonging to each string
		listIndex[i] = sampleHash.GIndex[list[i]]
	}
	return listIndex
}

func ByNames(inChan <-chan *Vcf, header *VcfHeader, list []string, writer *fileio.EasyWriter) {
	var listIndex []int16 = getListIndex(header, list)

	for record := range inChan {
		WriteVcf(writer, ReorderSampleColumns(record, listIndex))
	}
}

//FilterVcfPos will filter out records that appear as the same postion more than once, keeping the first one it encounters. In addition, if records contains Ns, those records will also be filtered out.
func FilterVcfPos(vcfs []*Vcf) []*Vcf {
	Sort(vcfs)
	var answer []*Vcf
	chrVcfMap := make(map[string][]*Vcf)

	var ref []dna.Base
	var alt [][]dna.Base
	var containsN bool
	for _, v := range vcfs {
		chrVcfMap[v.Chr] = append(chrVcfMap[v.Chr], v)
	}
	var i int
	var curr []*Vcf
	for key := range chrVcfMap {
		curr = chrVcfMap[key]
		encountered := make(map[int]bool)
		for i = 0; i < len(curr); i++ {
			if encountered[curr[i].Pos] == true {
				//do not add
			} else {
				encountered[curr[i].Pos] = true
				containsN = false
				ref = dna.StringToBases(curr[i].Ref)
				alt = GetAltBases(curr[i].Alt)
				if dna.CountBaseInterval(ref, dna.N, 0, len(ref)) != 0 {
					containsN = true
				}
				for j := 0; j < len(alt); j++ {
					if dna.CountBaseInterval(alt[j], dna.N, 0, len(alt[j])) != 0 {
						containsN = true
					}
				}
				if !containsN {
					answer = append(answer, curr[i])
				}
			}
		}
	}
	return answer
}
