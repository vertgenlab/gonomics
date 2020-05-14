package vcf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

func Filter(v *Vcf, chrom string, minPos int64, maxPos int64, ref string, alt string, minQual float64) bool {
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

func FilterQual(v *Vcf, minQual float64) bool {
	if v.Qual < minQual {
		return false
	}
	return true
}

func FilterAlt(v *Vcf, alt string) bool {
	if alt != "" && v.Alt != alt {
		return false
	}
	return true
}

func FilterRef(v *Vcf, ref string) bool {
	if ref != "" && v.Ref != ref {
		return false
	}
	return true
}

func FilterRange(v *Vcf, minPos int64, maxPos int64) bool {
	if v.Pos < minPos || v.Pos > maxPos {
		return false
	}
	return true
}

func FilterChrom(v *Vcf, chrom string) bool {
	if chrom != "" && v.Chr != chrom {
		return false
	}
	return true
}

func VcfOverlap(alpha *Vcf, beta *Vcf) bool {
	if common.MaxInt64(alpha.Pos-1, beta.Pos-1) < common.MinInt64(int64(len(dna.StringToBases(alpha.Ref))), int64(len(dna.StringToBases(beta.Ref)))) && strings.Compare(alpha.Chr, beta.Chr) == 0 {
		return true
	} else {
		return false
	}
}

func VcfBedOverlap(alpha *Vcf, beta *bed.Bed) bool {
	if common.MaxInt64(alpha.Pos-1, beta.ChromStart) < common.MinInt64(int64(len(dna.StringToBases(alpha.Ref))), beta.ChromEnd) && strings.Compare(alpha.Chr, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

func FilterAxtVcf(vcfs []*Vcf, fa []*fasta.Fasta) []*Vcf {
	split := VcfSplit(vcfs, fa)
	var answer []*Vcf
	var i, j int
	var ref []dna.Base
	var alt []dna.Base
	for i = 0; i < len(split); i++ {
		encountered := make(map[int64]bool)
		for j = 0; j < len(split[i]); j++ {
			if encountered[split[i][j].Pos] == true {
				//do not add
			} else {
				encountered[split[i][j].Pos] = true
				ref = dna.StringToBases(split[i][j].Ref)
				alt = dna.StringToBases(split[i][j].Alt)
				if dna.CountBaseInterval(ref, dna.N, 0, len(ref)) == 0 && dna.CountBaseInterval(alt, dna.N, 0, len(alt)) == 0 {
					answer = append(answer, split[i][j])
				}
			}
		}
	}
	Sort(answer)
	return answer
}

func FilterNs(vcfs []*Vcf) []*Vcf {
	var answer []*Vcf
	for i := 0; i < len(vcfs); i++ {
		if !strings.Contains(vcfs[i].Ref, "N") && !strings.Contains(vcfs[i].Alt, "N") {
			answer = append(answer, vcfs[i])
		}
	}
	return answer
}

func mergeSimilarVcf(a *Vcf, b *Vcf) *Vcf {
	mergeRecord := &Vcf{Chr: a.Chr, Pos: a.Pos, Id: a.Id, Ref: "", Alt: "", Qual: a.Qual, Filter: "Merged:SNP:INDEL", Info: a.Info, Format: "SVTYPE=SNP", Notes: a.Notes}
	if len(a.Ref) < len(b.Ref) {
		mergeRecord.Ref += b.Ref
	} else {
		mergeRecord.Ref += a.Ref
	}
	if len(a.Alt) < len(b.Alt) {
		mergeRecord.Alt += b.Alt
	} else {
		mergeRecord.Alt += a.Alt
	}
	return mergeRecord
}

func ApplyFilter(method string, gt []Haplotype, AaAa []int16, BBbb []int16) bool {
	switch {
	case strings.Contains(method, "AS"):
		if ASFilter(gt, AaAa, BBbb) {
			return true
		}
	case strings.Contains(method, "wrong"):
		if WrongFilter(gt, AaAa, BBbb) {
			return true
		}
	case strings.Contains(method, "medium"):
		if MediumHetFilter(gt, AaAa, BBbb) {
			return true
		}
	default:
		if ASFilter(gt, AaAa, BBbb) {
			return true
		}
	}
	return false
}

func ASFilter(gt []Haplotype, AaAa []int16, BBbb []int16) bool {
	if Heterozygous(AaAa, gt) && DifferentParents(gt[BBbb[0]], gt[BBbb[1]]) {
		return true
	} else {
		return false
	}
}

//get uninformative alleles: either parents are hets or the same allele
//use as bad training data for VQSR training
func WrongFilter(gt []Haplotype, AaAa []int16, BBbb []int16) bool {
	//gt := genotypeHelper(v)
	if Homozygous(AaAa, gt) || Heterozygous(BBbb, gt) {
		return true
	} else {
		return false
	}
}

func Homozygous(key []int16, gt []Haplotype) bool {
	for _, Aa := range key {
		if !IsHomozygous(gt[Aa]) {
			return false
		}
	}
	return true
}

func HetsOnly(gt []Haplotype, AaAa []int16) bool {
	//gt := genotypeHelper(v)
	if Heterozygous(AaAa, gt) {
		return true
	} else {
		return false
	}
}

//checks parents for homo alleles
//at least one het in other group plasses the filter.
func MediumHetFilter(samples []Haplotype, AaAa []int16, BBbb []int16) bool {
	//samples := genotypeHelper(v)
	if Homozygous(BBbb, samples) {
		for _, idx := range AaAa {
			if samples[idx].One != samples[idx].Two {
				return true
			}
		}
	}
	return false
}

func DifferentParents(parentOne Haplotype, parentTwo Haplotype) bool {
	if IsHomozygous(parentOne) && IsHomozygous(parentTwo) && (parentOne.One != parentTwo.One && parentOne.Two != parentTwo.Two) {
		return true
	} else {
		return false
	}
}

func WriteHeaderDev(file *fileio.EasyWriter, header *VcfHeader) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
	}
	common.ExitIfError(err)
}

func SimpleASFilter(v *Vcf, axb int16, a int16, b int16) bool {
	gt := GenotypeHelper(v)
	if LowQual(gt[axb]) || LowQual(gt[a]) || LowQual(gt[b]) {
		return false
	}
	if gt[axb].One != gt[axb].Two && (gt[a].One == gt[a].Two && (gt[a].One != gt[b].One || gt[a].Two != gt[b].Two)) && (gt[b].One == gt[b].Two) {
		return true
	} else {
		return false
	}
}
func ReadFilterList(filename string, dict map[string]int16) ([]string, []string) {
	file := fileio.Read(filename)
	if len(file) != 2 {
		log.Fatalf("Error: Sample sheet should be exactly 2 lines: hets one line 1, homo on line 2")
	}
	AaAa := strings.Split(file[0], ",")
	BBbb := strings.Split(file[1], ",")
	return AaAa, BBbb
	//return MapNameToIndex(dict, AaAa), MapNameToIndex(dict, BBbb)
}

func MapNameToIndex(dict map[string]int16, list []string) []int16 {
	var answer []int16 = make([]int16, len(list))
	for i := 0; i < len(list); i++ {
		answer[i] = dict[list[i]]
	}
	return answer
}

func LogEqual(alpha []*Vcf, beta []*Vcf) {
	if len(alpha) != len(beta) {
		log.Fatalf("len=%v and len=%v are not equal\n", len(alpha), len(beta))
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			log.Fatalf("%v and %v are not equal\n", alpha[i], beta[i])
		}
	}
}

func LowQual(hap Haplotype) bool {
	if hap.One < 0 || hap.Two < 0 {
		return true
	} else {
		return false
	}
}

func Unqualified(key []int16, gt []Haplotype) bool {
	for i := 0; i < len(key); i++ {
		if LowQual(gt[key[i]]) {
			return true
		}
	}
	return false
}

func IsHeterozygous(hap Haplotype) bool {
	if hap.One < 0 || hap.Two < 0 {
		return false
	}
	if hap.One != hap.Two {
		return true
	}
	if hap.One == hap.Two {
		return false
	}
	return false
}

func IsHomozygous(hap Haplotype) bool {
	if hap.One < 0 || hap.Two < 0 {
		return false
	}
	if hap.One == hap.Two {
		return true
	}
	if hap.One != hap.Two {
		return false
	}
	return false
}

func Heterozygous(key []int16, gt []Haplotype) bool {
	for _, Aa := range key {
		if !IsHeterozygous(gt[Aa]) {
			return false
		}
	}
	return true
}

//Both Homozygous, but different allele
func UniqueHomozygous(key []int16, AA []Haplotype) bool {
	hash := make(map[int16]bool)
	var ok bool
	for aa := 0; aa < len(key); aa++ {
		if IsHomozygous(AA[key[aa]]) {
			_, ok = hash[AA[key[aa]].One]
			if !ok {
				hash[AA[key[aa]].One] = true
			} else {
				return false
			}
		} else {
			return false
		}
	}
	return true
}

func StrictHetFilter(v *Vcf) bool {
	samples := GenotypeHelper(v)
	for _, gVcf := range samples {
		if gVcf.One == gVcf.Two {
			return false
		}
	}
	return true
}
