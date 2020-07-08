package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
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

func sameRecord(a *Vcf, b *Vcf) bool {
	if isEqual(a, b) {
		return true
	}
	if strings.Compare(a.Chr, b.Chr) == 0 && a.Pos == b.Pos {
		if strings.Compare(a.Ref, b.Ref) == 0 && strings.Compare(a.Alt, b.Alt) == 0 {
			return true
		}
	}
	return false
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

func ASFilter(v *Vcf, parentOne int16, parentTwo int16, F1 int16) bool {
	gt := GetAlleleGenotype(v)
	if isHomozygous(gt[parentOne]) && isHomozygous(gt[parentTwo]) && isHeterozygous(gt[F1]) {
		return true
	} else {
		return false
	}
}

func isHeterozygous(gt Sample) bool {
	if gt.AlleleOne < 0 || gt.AlleleTwo < 0 {
		return false
	}
	if gt.AlleleOne != gt.AlleleTwo {
		return true
	} else {
		return false
	}
}

func isHomozygous(gt Sample) bool {
	if gt.AlleleOne < 0 || gt.AlleleTwo < 0 {
		return false
	}
	if gt.AlleleOne == gt.AlleleTwo {
		return true
	} else {
		return false
	}
}

func Heterozygous(samples []int16, gt []Sample) bool {
	for _, Aa := range samples {
		if !isHeterozygous(gt[Aa]) {
			return false
		}
	}
	return true
}

func Homozygous(samples []int16, gt []Sample) bool {
	for _, Aa := range samples {
		if !isHomozygous(gt[Aa]) {
			return false
		}
	}
	return true
}
