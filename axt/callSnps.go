package axt

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"strconv"
	"strings"
)

func AxtToVcf(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	rCount := axtFile.RStart - 1 
	qCount := axtFile.QStart - 1
	for i := 0; i < len(axtFile.RSeq); i++ {

		var infoTag string
		if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), "-") != 0 && strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), "-") != 0 {
			rCount++
			qCount++
			//snp mismatch
			if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), dna.BaseToString(dna.ToUpper(axtFile.QSeq[i]))) != 0 {
				infoTag = "POS=" + strconv.FormatInt(qCount, 10)
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=SNP", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
				//fmt.Println(snps[i].RefSub, snps[i].QuerySub)
				answer = append(answer, curr)
			}
		}
		//insertion in VCF record
		if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), "-") == 0 {
			var altTmp string
			qCount++
			//var refTmp string
			altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.RSeq[j])), "-") == 0 {
					altTmp = altTmp + dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					qCount++
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=INS", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
		//deleteion vcf record
		if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), "-") == 0 {
			//var refTmp string
			var altTmp string
			tempRCount := 0
			altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if strings.Compare(dna.BaseToString(dna.ToUpper(axtFile.QSeq[j])), "-") == 0 {
					altTmp = altTmp + dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
					tempRCount++
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: altTmp, Alt: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=DEL", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
					rCount = rCount + int64(tempRCount)
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
	}
	return answer
}

func CallSnpsToVcf(axtList []*Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr []*vcf.Vcf
	for i := 0; i < len(axtList); i++ {
		curr = AxtToVcf(axtList[i])
		for _, j := range curr {
			answer = append(answer, j)
		}
	}
	return answer
}
//split vcf into slices to deal with different chromosomes
func VcfSplit(vcfRecord []*vcf.Vcf, fastaRecord []*fasta.Fasta) [][]*vcf.Vcf {
	var answer [][]*vcf.Vcf

	for i := 0; i < len(fastaRecord); i++ {
		var curr []*vcf.Vcf = nil
		for j :=0; j < len(vcfRecord); j++ {
			var pointer *vcf.Vcf
			if strings.Compare(fastaRecord[i].Name, vcfRecord[j].Chr ) == 0 {
				pointer = vcfRecord[j]
				curr = append(curr, pointer)
			}
		}
		vcf.Sort(curr)
		answer = append(answer, curr)
	}	
	return answer
}