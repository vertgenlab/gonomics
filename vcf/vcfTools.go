package vcf

import (
	"strings"
)

func FilterVcf(vcfs []*Vcf, snp bool, ins bool, del bool) []*Vcf {
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
