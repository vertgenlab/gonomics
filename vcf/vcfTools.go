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

func VcfToBed(v *Vcf) *bed.Bed {
	var b *bed.Bed = &bed.Bed{Chrom: v.Chr}
	if strings.Compare(v.Format, "SVTYPE=SNP") != 0 {
		b.ChromStart, b.ChromEnd = v.Pos-1, v.Pos-1
	}
	if strings.Compare(v.Format, "SVTYPE=INS") != 0 {
		b.ChromStart, b.ChromEnd = v.Pos, v.Pos
	}
	if strings.Compare(v.Format, "SVTYPE=DEL") != 0 {
		b.ChromStart = v.Pos
		b.ChromEnd = v.Pos + int64(len(dna.StringToBases(v.Ref))-1)
	}
	notes := strings.Join([]string{v.Format, v.Info, v.Notes}, "\t")
	b.Annotation = append(b.Annotation, notes)
	return b
}

func NoVcfOverlap(vcfs []*Vcf) {
	Sort(vcfs)
	var i, j int
	var alpha, beta *bed.Bed
	for i = 0; i < len(vcfs)-1; {
		alpha = VcfToBed(vcfs[i])
		beta = VcfToBed(vcfs[i+1])
		if !bed.Overlap(alpha, beta) {
			i++
		} else {
			for j = i + 1; j < len(vcfs)-1; j++ {
				vcfs[j] = vcfs[j+1]
			}
			vcfs = vcfs[:len(vcfs)-1]
		}
	}

	//for i := 0; i < len(vcfs)-1; i++ {
	//	if !strings.Contains(vcfs[i].Info, "SYTYPE=BND") {
	//		if !bed.Overlap(VcfToBed(vcfs[i]), VcfToBed(vcfs[i+1])) {
	//			answer = append(answer, vcfs[i])
	//		}
	//	}
	//}
	//return answer
}
