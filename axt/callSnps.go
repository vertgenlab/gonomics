package axt

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
	"fmt"
)

func AxtToVcf(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	rCount := axtFile.RStart - 1
	qCount := axtFile.QStart - 1
	for i := 0; i < len(axtFile.RSeq); i++ {
		if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
			rCount++
			qCount++
			//snp mismatch
			if dna.ToUpper(axtFile.RSeq[i]) != dna.ToUpper(axtFile.QSeq[i]) {
				if dna.ToUpper(axtFile.RSeq[i]) != dna.N || dna.ToUpper(axtFile.QSeq[i]) != dna.N {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 0, Filter: "PASS", Info: fmt.Sprint(qCount), Format: "SVTYPE=SNP", Notes: "noNs"}
					answer = append(answer, curr)
				} else if dna.ToUpper(axtFile.RSeq[i]) == dna.N && dna.ToUpper(axtFile.QSeq[i]) != dna.N {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]), Alt: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Qual: 0, Filter: "PASS", Info: fmt.Sprint(qCount), Format: "SVTYPE=INS", Notes: "Closed gap in reference"}
					for j := i; j < len(axtFile.RSeq); j++ {
						if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
							rCount++
							qCount++
							if dna.ToUpper(axtFile.RSeq[j]) == dna.N && dna.ToUpper(axtFile.QSeq[j]) != dna.N {
								curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
							} else {
								answer = append(answer, curr)
								i = j - 1
								break
							}
						}
					}
				} else if dna.ToUpper(axtFile.QSeq[i]) == dna.N && dna.ToUpper(axtFile.RSeq[i]) != dna.N {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]), Alt: dna.BaseToString(axtFile.RSeq[i-1]), Qual: 0, Filter: "PASS", Info: fmt.Sprint(qCount), Format: "SVTYPE=DEL", Notes: "Closed gap query"}
					for j := i; j < len(axtFile.RSeq); j++ {
						if dna.ToUpper(axtFile.RSeq[j]) != dna.ToUpper(axtFile.QSeq[j]) {
							rCount++
							qCount++
							if dna.ToUpper(axtFile.QSeq[j]) == dna.N && dna.ToUpper(axtFile.RSeq[j]) != dna.N {
								curr.Ref += dna.BaseToString(axtFile.RSeq[j])
							} else {
								answer = append(answer, curr)
								i = j - 1
								break
							}
						}

					}
				}
			}
		}
		//insertion in VCF record
		if axtFile.RSeq[i] == dna.Gap {
			var altTmp string
			qCount++
			//var refTmp string
			altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if axtFile.RSeq[j] == dna.Gap {
					if dna.ToUpper(axtFile.QSeq[j]) != dna.N {
						altTmp = altTmp + dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					}
					qCount++
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: fmt.Sprint(qCount), Format: "SVTYPE=INS", Notes: "noNs"}
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
		//deleteion vcf record
		if axtFile.QSeq[i] == dna.Gap {
			//var refTmp string
			var altTmp string
			tempRCount := 0
			altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if axtFile.QSeq[j] == dna.Gap {
					if dna.ToUpper(axtFile.RSeq[j]) != dna.N {
						altTmp = altTmp + dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
					}
					tempRCount++
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: altTmp, Alt: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Qual: 0, Filter: "PASS", Info: fmt.Sprint(qCount), Format: "SVTYPE=DEL", Notes: "noNs"}
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

/*
func AxtToVcf(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	rCount := axtFile.RStart - 1
	qCount := axtFile.QStart - 1
	var i, j int
	dna.AllToUpper(axtFile.RSeq)
	dna.AllToUpper(axtFile.QSeq)
	for i = 0; i < len(axtFile.RSeq); i++ {
		if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
			rCount++
			qCount++
			//snp mismatch
			if (axtFile.RSeq[i]!= axtFile.QSeq[i]) && (axtFile.RSeq[i] != dna.N) && (axtFile.QSeq[j] != dna.N) {

				//if strings.Compare(dna.BaseToString(axtFile.RSeq[i]), dna.BaseToString(axtFile.QSeq[i])) != 0 {
				//infoTag = "POS=" + strconv.FormatInt(qCount, 10)
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i]), Alt: dna.BaseToString(axtFile.QSeq[i]), Qual: 0, Filter: "", Format: "SVTYPE=SNP", Notes: ""}
				//fmt.Println(snps[i].RefSub, snps[i].QuerySub)
				answer = append(answer, curr)
			}
			if (axtFile.RSeq[i] == dna.N) && (axtFile.QSeq[i] != dna.N) {
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]), Alt: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))+dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 0, Filter: "", Info: "", Format: "SVTYPE=INS", Notes: "Closed gap in reference"}
				for j = 1; j < len(axtFile.RSeq); j++ {
					if (axtFile.QSeq[i+j] != dna.Gap) && (axtFile.QSeq[i+j] != dna.N) && (axtFile.RSeq[i+j] != dna.Gap) && (axtFile.RSeq[i+j] == dna.N) {
						curr.Alt += dna.BaseToString(axtFile.QSeq[j])
					} else {
						answer = append(answer, curr)
						rCount = rCount + int64(j)
						qCount = qCount + int64(j)
						i = j - 1
						break
					}
				}
			}
			if (axtFile.RSeq[i] != dna.N) && (axtFile.QSeq[i] == dna.N) {
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]) + dna.BaseToString(axtFile.RSeq[i]), Alt: dna.BaseToString(axtFile.RSeq[i-1]), Qual: 0, Filter: "", Info: "", Format: "SVTYPE=DEL", Notes: "Closed gap query"}
				for j = 1; j < len(axtFile.RSeq); j++ {
					if axtFile.RSeq[i+j] != dna.Gap && axtFile.RSeq[i+j] != dna.N && (axtFile.QSeq[i+j] != dna.Gap) && (axtFile.QSeq[i+j] == dna.N) {
						curr.Ref += dna.BaseToString(axtFile.RSeq[i+j])
					} else {
						answer = append(answer, curr)
						rCount = rCount + int64(j)
						qCount = qCount + int64(j)
						i = j -1
						break
					}

				}
			}

		}
		//insertion in VCF record
		if axtFile.RSeq[i] == dna.Gap && (axtFile.QSeq[i] != dna.Gap || axtFile.QSeq[i] == dna.N) {
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]), Alt: dna.BaseToString(axtFile.RSeq[i-1])+ dna.BaseToString(axtFile.QSeq[i]), Qual: 0, Filter: "", Info: "", Format: "SVTYPE=INS", Notes: ""}
			for j = i; j < len(axtFile.RSeq); j++ {
				if axtFile.RSeq[j] == dna.Gap && (axtFile.QSeq[j] != dna.Gap || axtFile.QSeq[j] == dna.N) {
					qCount++
					curr.Alt += dna.BaseToString(axtFile.QSeq[i+j])
				} else {
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
		//deleteion vcf record
		if axtFile.QSeq[i] == dna.Gap {
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]) + dna.BaseToString(axtFile.RSeq[i]), Alt: dna.BaseToString(axtFile.RSeq[i-1]), Qual: 0, Filter: "", Info: "", Format: "SVTYPE=DEL", Notes: ""}
			for j = i; j < len(axtFile.RSeq); j++ {
				if axtFile.QSeq[j] == dna.Gap && (axtFile.RSeq[j] != dna.Gap || axtFile.RSeq[j] == dna.N) {
					rCount++
					curr.Ref += dna.BaseToString(axtFile.RSeq[j])
				} else {
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
	}
	return answer
}*/
