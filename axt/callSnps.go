package axt

import (
	"fmt"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
)

// ToVcfFile takes a
func ToVcfFile(filename string, axtList []Axt, fa []fasta.Fasta) {
	var records []vcf.Vcf
	for i := range axtList {
		records = append(records, ToVcf(axtList[i])...)
	}
	vcf.Write(filename, records)
}

func info(input Axt) string {
	var text string = ""
	text = fmt.Sprintf("%s;%d;%d;%s;%d;%d;%t;%d", input.RName, input.RStart, input.REnd, input.QName, input.QStart, input.QEnd, input.QStrandPos, input.Score)
	return text
}

func ToVcf(axtFile Axt) []vcf.Vcf {
	var answer []vcf.Vcf
	var curr vcf.Vcf
	var rCount int = axtFile.RStart - 1
	qCount := axtFile.QStart - 1
	for i := 0; i < len(axtFile.RSeq); i++ {
		if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
			rCount++
			qCount++
			//snp mismatch
			if dna.ToUpper(axtFile.RSeq[i]) != dna.ToUpper(axtFile.QSeq[i]) {
				curr = vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: []string{dna.BaseToString(dna.ToUpper(axtFile.QSeq[i]))}, Qual: 30, Filter: "PASS", Info: fmt.Sprintf("query=%d;SVTYPE=SNP;%s", qCount, info(axtFile))}
				answer = append(answer, curr)
			}
		}
		//insertion in VCF record
		if axtFile.RSeq[i] == dna.Gap {

			qCount++
			curr = vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: []string{dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1]))}, Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d;SVTYPE=SNP;%s", qCount, info(axtFile))}

			for j := i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.RSeq[j]) == dna.Gap {
					curr.Alt[0] += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					qCount++
				} else {
					if len(answer) == 0 {
						answer = append(answer, curr)
					} else if answer[len(answer)-1].Pos == curr.Pos && strings.Compare(answer[len(answer)-1].Info, "SVTYPE=SNP") == 0 {
						curr.Info = "SVTYPE=SNP;INS"
						answer[len(answer)-1] = curr
					} else {
						answer = append(answer, curr)
					}
					//curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=INS", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
		//deleteion vcf record
		if axtFile.QSeq[i] == dna.Gap {
			tempRCount := 0
			curr = vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: []string{dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1]))}, Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d;SVTYPE=DEL;%s", qCount, info(axtFile))}
			//altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.QSeq[j]) == dna.Gap {
					curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
					//altTmp = altTmp +
					tempRCount++
				} else {
					if len(answer) == 0 {
						answer = append(answer, curr)
					} else if answer[len(answer)-1].Pos == curr.Pos && strings.Compare(answer[len(answer)-1].Info, "SVTYPE=SNP") == 0 {
						curr.Info = "SVTYPE=SNP;DEL"
						answer[len(answer)-1] = curr
					} else {
						answer = append(answer, curr)
					}
					rCount = rCount + tempRCount
					i = j - 1
					break
				}
			}
		}
	}
	//log.Printf("\nFound %d differences in this block...\n%s\n", len(answer), AxtInfo(axtFile))
	return answer
}
