package vcf

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"strings"
	

)

func SelectVcf(vcfs []*Vcf, snp bool, ins bool, del bool) []*Vcf {
	var answer []*Vcf
	for i := 0; i < len(vcfs); i++ {
		if !snp {
			if strings.Compare(vcfs[i].Format, "SVTYPE=SNP") != 0 {
				answer = append(answer, vcfs[i])
			}
		} else if !ins {
			if strings.Compare(vcfs[i].Format, "SVTYPE=INS") != 0 {
				answer = append(answer, vcfs[i])
			}
		} else if !del {
			if strings.Compare(vcfs[i].Format, "SVTYPE=DEL") != 0 {
				answer = append(answer, vcfs[i])
			}
		}
	}
	return answer
}

func VcfToBed(vcfs []*Vcf) []*bed.Bed {
	var answer []*bed.Bed = make([]*bed.Bed, len(vcfs))
	for i := 0; i < len(vcfs); i++ {
		answer[i] = vcfLineToBed(vcfs[i])
	}
	return answer
}

func vcfLineToBed(v *Vcf) *bed.Bed {
	var b *bed.Bed = &bed.Bed{Chrom: v.Chr, Name: v.Id}
	if  Snp(v) {
		b.ChromStart, b.ChromEnd = v.Pos-1, v.Pos-1
	}
	if  Ins(v) {
		b.ChromStart, b.ChromEnd = v.Pos, v.Pos
	}
	if  Del(v) {
		b.ChromStart = v.Pos
		b.ChromEnd = v.Pos + int64(len(dna.StringToBases(v.Ref))-1)
	}
	return b
}

func NoVcfOverlap(vcfs []*Vcf) {
	Sort(vcfs)
	var i, j int
	var alpha, beta *bed.Bed
	for i = 0; i < len(vcfs)-1; {
		alpha = vcfLineToBed(vcfs[i])
		beta = vcfLineToBed(vcfs[i+1])
		if !bed.Overlap(alpha, beta) {
			i++
		} else {
			for j = i + 1; j < len(vcfs)-1; j++ {
				vcfs[j] = vcfs[j+1]
			}
			vcfs = vcfs[:len(vcfs)-1]
		}
	}

}

func Snp(v *Vcf) bool {
	if strings.Compare(v.Format, "SVTYPE=SNP") == 0 {
		return true
	}
	if strings.Compare(v.Info, "SVTYPE=SNP") == 0 {
		return true
	}
	return false
}

func Ins(v *Vcf) bool {
	if strings.Compare(v.Format, "SVTYPE=INS") == 0 {
		return true
	}
	if strings.Contains(v.Info, "SVTYPE=INS") {
		return true
	}
	return false
}

func Del(v *Vcf) bool {
	var truth bool = false
	if strings.Compare(v.Format, "SVTYPE=DEL") == 0 {
		return true
	}
	if strings.Contains(v.Info, "SVTYPE=DEL") {
		return true
	}
	return truth
}

func CopyVcfPointer(v *Vcf) *Vcf {
	answer := &Vcf{Chr: v.Chr, Pos: v.Pos, Id: v.Id, Ref: v.Ref, Alt: v.Alt, Qual: v.Qual, Filter: v.Filter, Info: v.Info, Format: v.Format, Notes: v.Notes}
	return answer
}
/*
func vcfMergeAllelles(a *Vcf, b *Vcf, aXb *Vcf) {
	if strings.Compare(a.Chr, b.Chr) == 0 && a.Pos == b.Pos {
		log.Printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", a.Chr, a.Pos, a.Alt, b.Alt, aXb.Alt, a.Notes, b.Notes, aXb.Notes)
	}
}

func VcfToAlleleMap(a []*Vcf, b []*Vcf, aXb []*Vcf) {
	for i := 0; i < common.Min(common.Min(len(a), len(b)), len(aXb)); i++ {
		vcfMergeAllelles(a[i], b[i], aXb[i])
	}
}*/


/*
type VcfAllele struct {
	Chr string
	Pos int64
	Hybrid *Info
	A1 *Info
	A2 *Info
}

type Info struct {
	Ref string
	Alt string
	Notes string
}

func vcfToInfoHelp(v *Vcf) *Info {
	answer := Info{Ref: v.Chr, Alt: v.Alt, Notes: v.Notes}
	return &answer
}*/
/*
func MergeVcfAlleles(aXb []*Vcf, a []*Vcf, b []*Vcf) []*VcfAllele {
	var answer []*VcfAllele
	alleleOne := VcfToMap(a)
	alleleTwo := VcfToMap(b)
	var lookUp string
	for _,v := range aXb {
		lookUp = chromPosToString(v.Chr, v.Pos)
		_, okOne := alleleOne[lookUp]
		_, okTwo := alleleTwo[lookUp]
		if  {
			curr := &VcfAllele{Chr: v.Chr, Pos: v.Pos,  Hybrid: vcfToInfoHelp(v), A1: vcfToInfoHelp(alleleOne[lookUp]), A2: vcfToInfoHelp(alleleTwo[lookUp])}
			answer = append(answer, curr)
		}
	}
	return answer
}*/
/*
func MergeVcf(a []*Vcf, b []*Vcf) []*Vcf {
	var answer []*Vcf
	var i int
	var lookUp *Vcf
	mapToVcf := VcfToMap(b)
	for i = 0; i < len(a); i++ {
		lookUp = mapToVcf[chromPosToString(a[i].Chr, a[i].Pos)]
		if lookUp != nil {
			if strings.Compare(a[i].Chr, lookUp.Chr) == 0 && a[i].Pos == lookUp.Pos {
				curr := CopyVcfPointer(a[i])
				curr.Notes = fmt.Sprintf("%s\t%s", a[i].Notes, lookUp.Notes)
				answer = append(answer, curr)
			}
		}

	}
	return answer
}*/