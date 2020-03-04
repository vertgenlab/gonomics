package axt

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
)

/*
func CallSnpsToVcf(axtList []*Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr []*vcf.Vcf
	for i := 0; i < len(axtList); i++ {
		curr = HaplotypeBlocks(axtList[i])
		vcf.Sort(curr)
		answer = append(answer, curr...)
	}
	return answer
}*/

func AxtVcfToFile(filename string, axtList []*Axt, fa []*fasta.Fasta) {

	var records []*vcf.Vcf
	for i := 0; i < len(axtList); i++ {
		records = append(records, AxtToVcf(axtList[i])...)
	}
	records = vcf.FilterAxtVcf(records, fa)

	//vcfs := vcf.Read(filename)
	//sorted := fileio.MustCreate(filename + ".sorted.vcf")
	//vcf.WriteHeader(records, head)
	vcf.Write(filename, records)
}

func AxtGapsVcfToFile(filename string, axtList []*Axt, fa []*fasta.Fasta) {
	ref := fasta.FastaMap(fa)
	var records []*vcf.Vcf
	var refIndex int64 = 0
	var lastChr string = axtList[0].RName
	var gap *vcf.Vcf
	var refSeq []dna.Base
	for i := 0; i < len(axtList); i++ {
		if axtList[i].RStart - refIndex > 1 && strings.Compare(lastChr, axtList[i].RName) == 0 {
			refSeq = ref[axtList[i].RName][refIndex:axtList[i].RStart]
			dna.AllToUpper(refSeq)
			gap = &vcf.Vcf{Chr: axtList[i].RName, Pos: refIndex, Id: axtList[i].QName, Ref: dna.BasesToString(refSeq), Alt: dna.BaseToString(dna.ToUpper(ref[axtList[i].RName][refIndex])), Qual: 24, Filter: "PASS", Info: ".", Format: "SVTYPE=DEL", Notes: AxtInfo(axtList[i])}
			records = append(records, gap)
		}
		records = append(records, AxtToVcf(axtList[i])...)
		refIndex = axtList[i].REnd-1
		lastChr = axtList[i].RName
	}
	records = vcf.FilterAxtVcf(records, fa)

	//vcfs := vcf.Read(filename)
	//sorted := fileio.MustCreate(filename + ".sorted.vcf")
	//vcf.WriteHeader(records, head)
	vcf.Write(filename, records)
}

/*
func FindBiggerMutations(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	//axtFile.RStart
	//axtFile.QStart
	var target int64 = axtFile.RStart
	var tGaps int64 = 0
	var query int64 = axtFile.QStart
	var qGaps int64 = 0
	var i, j int
	for i = 0; i < len(axtFile.RSeq); i++ {
		if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
			//snp mismatch
			if dna.ToUpper(axtFile.RSeq[i]) != dna.ToUpper(axtFile.QSeq[i]) {
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: target+int64(i)-tGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 30, Filter: "PASS", Info: fmt.Sprintf("query=%d", query+int64(i)-qGaps), Format: "SVTYPE=SNP", Notes: AxtInfo(axtFile)}
				//loop ahead to check for snps next to each other
				for j = i + 1; j < len(axtFile.RSeq); j++ {
					if axtFile.RSeq[j] != dna.Gap && axtFile.QSeq[j] != dna.Gap && dna.ToUpper(axtFile.RSeq[i]) != dna.ToUpper(axtFile.QSeq[i]) {
						curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
						curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					} else {
						answer = append(answer, curr)
						i = j - 1
						break
					}
				}
			}
		}
		//insertion in VCF record
		if axtFile.RSeq[i] == dna.Gap {
			tGaps++
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: target+int64(i)-tGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", query+int64(i)-qGaps), Format: "SVTYPE=INS", Notes: AxtInfo(axtFile)}
			curr.Alt +=dna.BaseToString(dna.ToUpper(axtFile.QSeq[i]))
			for j = i + 1; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.RSeq[j]) == dna.Gap && axtFile.QSeq[j] != dna.Gap {
					curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
				} else {
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
		//deleteion vcf record
		if axtFile.QSeq[i] == dna.Gap  {
			qGaps++
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: target+int64(i-1)-tGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", query+int64(i)-qGaps), Format: "SVTYPE=DEL", Notes: AxtInfo(axtFile)}
			for j = i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.QSeq[j]) == dna.Gap {
					curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))

				} else {
					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
	}
	log.Printf("\nFound %d differences in this block...\n%s\n", len(answer), AxtInfo(axtFile))
	return answer
}*/
/*
func FindSmallMutations(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	dna.AllToUpper(axtFile.RSeq)
	dna.AllToUpper(axtFile.QSeq)
	var targetStart int64 = axtFile.RStart
	var targetGaps int64 = 0
	var queryStart int64 = axtFile.QStart
	var queryGaps int64 = 0
	var i, j int
	for i = 0; i < len(axtFile.RSeq); i++ {
		if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
			//snp mismatch
			if i < len(axtFile.RSeq)-1 {
				if axtFile.RSeq[i+1] == dna.Gap || axtFile.QSeq[i+1] == dna.Gap {
					continue
				}
			}
			if axtFile.RSeq[i] != axtFile.QSeq[i] && axtFile.RSeq[i] != dna.N {
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: axtFile.RStart + int64(i) - targetGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 30, Filter: "PASS", Info: fmt.Sprintf("query=%d", queryStart+int64(i)-queryGaps), Format: "SVTYPE=SNP", Notes: AxtInfo(axtFile)}
				answer = append(answer, curr)
			}
		} else if axtFile.RSeq[i] == dna.Gap {
			//insertion in VCF record
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: targetStart + int64(i) - targetGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", queryStart+int64(i)-queryGaps), Format: "SVTYPE=INS", Notes: AxtInfo(axtFile)}
			for j = i; j < len(axtFile.RSeq); j++ {
				if axtFile.RSeq[j] == dna.Gap {
					curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					queryGaps++
				} else {

					answer = append(answer, curr)

					i = j - 1
					break
				}
			}
		} else if axtFile.QSeq[i] == dna.Gap {
			//deleteion vcf record
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: targetStart + int64(i-1) - targetGaps, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", queryStart+int64(i)-queryGaps), Format: "SVTYPE=DEL", Notes: AxtInfo(axtFile)}
			for j = i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.QSeq[j]) == dna.Gap {
					curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
					targetGaps++
				} else {

					answer = append(answer, curr)
					i = j - 1
					break
				}
			}
		}
	}
	log.Printf("\nFound %d differences in this block...\n%s\n", len(answer), AxtInfo(axtFile))
	return answer
}

func HaplotypeBlocks(axtFile *Axt) []*vcf.Vcf {
	var curr *vcf.Vcf
	var answer []*vcf.Vcf
	rCount := axtFile.RStart - 1
	qCount := axtFile.QStart - 1
	var i, j int
	for i = 0; i < len(axtFile.RSeq); i++ {
		if axtFile.QSeq[i] == dna.N || axtFile.RSeq[i] == dna.N {
			if axtFile.RSeq[i] == dna.N {
				rCount++
			}
			if axtFile.QSeq[i] == dna.N {
				qCount++
			}
		} else {
			if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
				rCount++
				qCount++
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: "", Alt: "", Qual: 24, Filter: "PASS", Info: AxtInfo(axtFile), Format: "SVTYPE=SNP", Notes: "RABStoGasAcu1-4"}
				for j = i; j < len(axtFile.RSeq); j++ {
					if axtFile.RSeq[i] != dna.Gap && axtFile.QSeq[i] != dna.Gap {
						if (axtFile.RSeq[i] != dna.N) && (axtFile.QSeq[i] != dna.N) && (dna.ToUpper(axtFile.RSeq[i]) != dna.ToUpper(axtFile.QSeq[i])) {
							rCount++
							qCount++
							curr.Notes += fmt.Sprintf("QueryPos=%d", qCount)
							curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
							curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
						} else {
							//vcf.WriteVcf(file, curr)
							answer = append(answer, curr)
							i = j
							break
						}
					}
				}
			} else if axtFile.RSeq[i] == dna.Gap { //insertion in VCF record
				qCount++
				if dna.ToUpper(axtFile.RSeq[i-1]) == dna.ToUpper(axtFile.QSeq[i-1]) && dna.ToUpper(axtFile.RSeq[i-1]) != dna.N {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(axtFile.QSeq[i-1]), Qual: 24, Filter: "PASS", Info: fmt.Sprintf(";queryPos=%d", qCount), Format: "SVTYPE=INS", Notes: AxtInfo(axtFile)}
					for j = i; j < len(axtFile.RSeq); j++ {
						if axtFile.RSeq[i] == dna.Gap && axtFile.QSeq[i] != dna.Gap && axtFile.QSeq[j] != dna.N {
							curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
							qCount++
						} else {
							//vcf.WriteVcf(file, curr)
							answer = append(answer, curr)
							i = j
							break
						}
					}
				}
			} else if axtFile.QSeq[i] == dna.Gap { //deleteion vcf record
				if dna.ToUpper(axtFile.RSeq[i-1]) == dna.ToUpper(axtFile.QSeq[i-1]) && dna.ToUpper(axtFile.RSeq[i-1]) != dna.N {
					curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("QueryPos=%d", qCount), Format: "SVTYPE=DEL", Notes: AxtInfo(axtFile)}
					for j = i; j < len(axtFile.RSeq); j++ {
						if axtFile.RSeq[j] != dna.Gap && axtFile.QSeq[j] == dna.Gap && axtFile.RSeq[j] != dna.N {
							rCount++
							curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
						} else {
							//vcf.WriteVcf(file, curr)
							answer = append(answer, curr)
							i = j
							break
						}
					}
				}
			}
		}
	}
	//log.Printf("\nFound %d differences in this block...\n", AxtInfo(axtFile), len(answer))
	return answer
}*/

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
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 30, Filter: "PASS", Info: fmt.Sprintf("query=%d", qCount), Format: "SVTYPE=SNP", Notes: AxtInfo(axtFile)}
				answer = append(answer, curr)
			}
		}
		//insertion in VCF record
		if axtFile.RSeq[i] == dna.Gap {

			qCount++
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", qCount), Format: "SVTYPE=INS", Notes: AxtInfo(axtFile)}

			for j := i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.RSeq[j]) == dna.Gap {
					curr.Alt += dna.BaseToString(dna.ToUpper(axtFile.QSeq[j]))
					qCount++
				} else {
					if len(answer) == 0 {
						answer = append(answer, curr)
					} else if answer[len(answer)-1].Pos == curr.Pos && strings.Compare(answer[len(answer)-1].Format, "SVTYPE=SNP") == 0 {
						curr.Format = "SVTYPE=SNP;INS"
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
			curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])), Alt: dna.BaseToString(dna.ToUpper(axtFile.QSeq[i-1])), Qual: 24, Filter: "PASS", Info: fmt.Sprintf("query=%d", qCount), Format: "SVTYPE=DEL", Notes: AxtInfo(axtFile)}
			//altTmp = dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1]))
			for j := i; j < len(axtFile.RSeq); j++ {
				if dna.ToUpper(axtFile.QSeq[j]) == dna.Gap {
					curr.Ref += dna.BaseToString(dna.ToUpper(axtFile.RSeq[j]))
					//altTmp = altTmp +
					tempRCount++
				} else {
					if len(answer) == 0 {
						answer = append(answer, curr)
					} else if answer[len(answer)-1].Pos == curr.Pos && strings.Compare(answer[len(answer)-1].Format, "SVTYPE=SNP") == 0 {
						curr.Format = "SVTYPE=SNP;DEL"
						answer[len(answer)-1] = curr
					} else {
						answer = append(answer, curr)
					}
					rCount = rCount + int64(tempRCount)
					i = j - 1
					break
				}
			}
		}
	}
	log.Printf("\nFound %d differences in this block...\n%s\n", len(answer), AxtInfo(axtFile))
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
			if (axtFile.RSeq[i] != axtFile.QSeq[i]) && (axtFile.RSeq[i] != dna.N) && (axtFile.QSeq[j] != dna.N) {

				//if strings.Compare(dna.BaseToString(axtFile.RSeq[i]), dna.BaseToString(axtFile.QSeq[i])) != 0 {
				//infoTag = "POS=" + strconv.FormatInt(qCount, 10)
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i]), Alt: dna.BaseToString(axtFile.QSeq[i]), Qual: 0, Filter: "", Format: "SVTYPE=SNP", Notes: ""}
				//fmt.Println(snps[i].RefSub, snps[i].QuerySub)
				answer = append(answer, curr)
			}
			if (axtFile.RSeq[i] == dna.N) && (axtFile.QSeq[i] != dna.N) {
				curr = &vcf.Vcf{Chr: axtFile.RName, Pos: rCount, Id: axtFile.QName, Ref: dna.BaseToString(axtFile.RSeq[i-1]), Alt: dna.BaseToString(dna.ToUpper(axtFile.RSeq[i-1])) + dna.BaseToString(dna.ToUpper(axtFile.QSeq[i])), Qual: 0, Filter: "", Info: "", Format: "SVTYPE=INS", Notes: "Closed gap in reference"}
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
						i = j - 1
						break
					}

				}
			}

		}

	}
	return answer
}*/

/*
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
}*/
