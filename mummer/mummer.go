package mummer

import (
	"os"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vergenlab/gonomics/vcf"
	)

type Snp struct {
	RefName   string
	RefPos    int64
	RefSub    string
	QueryName string
	QueryPos  int64
	QuerySub  string
}

func MummerToVcf(mummer []*Snp, reference map[string][]dna.Base) []*Vcf {
	var answer []*Vcf
	var curr *Vcf
	var currentSeq []dna.Base
	for i := 0; i < len(mummer); i++ {
		switch {
		//SNP
		case SnpTruth(mummer[i]) == 0:
			curr = &Vcf{Chr: mummer[i].RefName, Pos: mummer[i].RefPos, Id: ".", Ref: mummer[i].RefSub, Alt: mummer[i].QuerySub, Qual: 0, Filter: "PASS", Info: ".", Format: "IREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=0;QR=0;RO=0;RPL=0;RPP=7.35324;RPPR=0;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=0;SRP=0;SRR=0;TYPE=snp", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
			//fmt.Println(mummer[i].RefSub, mummer[i].QuerySub)
			answer = append(answer, curr)
		//logic for insertion relative to the reference
		case strings.Compare(mummer[i].RefSub, ".") == 0:
			var altTmp string
			currentSeq = reference[mummer[i].RefName]
			altTmp = dna.BaseToString(dna.ToUpper(currentSeq[mummer[i].RefPos-1])) + mummer[i].QuerySub
			curr = &Vcf{Chr: mummer[i].RefName, Pos: mummer[i].RefPos, Id: ".", Ref: dna.BaseToString(dna.ToUpper(currentSeq[mummer[i].RefPos-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: ".", Format: "IREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=0;QR=0;RO=0;RPL=0;RPP=7.35324;RPPR=0;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=0;SRP=0;SRR=0;TYPE=snp", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
			answer = append(answer, curr)
		//logic for deletions
		case strings.Compare(mummer[i].RefSub, ".") == 0:
			var refTmp string
			currentSeq = reference[mummer[i].RefName]
			refTmp = dna.BaseToString(dna.ToUpper(currentSeq[mummer[i].RefPos-1])) + dna.BaseToString(dna.ToUpper(currentSeq[mummer[i].RefPos-1]))
			//var refTmp string
			curr = &Vcf{Chr: mummer[i].RefName, Pos: mummer[i].RefPos, Id: ".", Ref: refTmp, Alt: dna.BaseToString(dna.ToUpper(currentSeq[mummer[i].RefPos-1])), Qual: 0, Filter: "PASS", Info: ".", Format: "MQMR=0;NS=1;NUMALT=1;ODDS=4.15888;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=0;QR=0;RO=0;RPL=2;RPP=7.35324;RPPR=0;RPR=0;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=0;SRP=0;SRR=0;TYPE=del", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
			answer = append(answer, curr)
		}
	}
	return answer
}

func PrintSnp(input []*Snp) {
	for i := range input {
		fmt.Printf("%s\t%s\t%v\t%v\t%s\t%s\n", input[i].RefName, input[i].QueryName, input[i].RefPos, input[i].QueryPos, input[i].RefSub, input[i].QuerySub)
	}
}

func MummerSNP(mummerFile string) ([]*Snp, error) {
	var answer []*Snp
	var curr *Snp
	var line string
	file, err := os.Open(mummerFile)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, "\t")
		//fmt.Println("there is data here")
		switch {
		case len(data) == 1:
			//these lines are sequences, and we are not recording them
			//fmt.Println("found sequences")
		case len(line) == 0:
			//blank line
			//fmt.Println("found blank")
		case len(data) == 12:
			//fmt.Println("found header line")

			pRef, _ := strconv.ParseInt(data[0], 10, 64)
			pQuery, _ := strconv.ParseInt(data[3], 10, 64)
			curr = &Snp{RefName: data[10], RefPos: pRef, RefSub: data[1], QueryName: data[11], QueryPos: pQuery, QuerySub: data[2]}
			answer = append(answer, curr)

		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer, nil
}


func SnpTruth(s *Snp) int {
	if strings.Compare(s.RefSub, s.QuerySub) != 0 && strings.Compare(s.RefSub, ".") != 0 && strings.Compare(s.QuerySub, ".") != 0 {
		return 0
	} else {
		return -1
	}
}