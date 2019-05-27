package axt

import (
	"bufio"
	//"flag"
	"fmt"
	"sort"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"io/ioutil"
	"log"
	"os"
	"strconv"
	"strings"
)

type Axt struct {
	//Next *Axt
	AlignNumber   int64
	RefName       string
	RefStart      int64
	RefEnd        int64
	AligningName  string
	AligningStart int64
	AligningEnd   int64
	Strand        string
	BlastzScore   int64
	RefSeq        string
	QuerySeq      string
}

func ReadIn(filename string) ([]*Axt, error) {
	var answer []*Axt
	var curr *Axt
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	var prefix bool
	rline, prefix, err2 = reader.ReadLine()
	for ; err2 != io.EOF; rline, prefix, err2 = reader.ReadLine() {
		line = string(rline[:])
		//header data:
		for strings.HasPrefix(line, "#") {
			rline, prefix, err2 = reader.ReadLine()
			line = string(rline[:])
		}
		data := strings.Split(line, " ")
		var an int64
		var rs int64
		var re int64
		var as int64
		var ae int64
		var bs int64
		var seq1 string
		var seq2 string

		if len(data) != 9 {
			fmt.Println("something went wrong")
		}
		an, _ = strconv.ParseInt(data[0], 10, 64)
		rs, _ = strconv.ParseInt(data[2], 10, 64)
		re, _ = strconv.ParseInt(data[3], 10, 64)
		as, _ = strconv.ParseInt(data[5], 10, 64)
		ae, _ = strconv.ParseInt(data[6], 10, 64)
		bs, _ = strconv.ParseInt(data[8], 10, 64)

		rline, prefix, err2 = reader.ReadLine()
		seq1 = string(rline[:])
		for prefix != false {
			rline, prefix, err2 = reader.ReadLine()
			if err2 != nil {
				log.Fatal(err2)
			}
			seq1 += string(rline[:])
		}

		rline, prefix, err2 = reader.ReadLine()
		seq2 = string(rline[:])
		for prefix != false {
			rline, prefix, err2 = reader.ReadLine()
			if err2 != nil {
				log.Fatal(err2)
			}
			line = string(rline[:])
			seq2 += line
		}
		if len(seq1) != len(seq2) {
			fmt.Println("Sequences are not the same length")
			log.Fatal(err2)
		}

		rline, prefix, err2 = reader.ReadLine()
		blank := string(rline[:])
		if len(blank) != 0 {
			fmt.Println("not a blank line")
		}
		curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: as, AligningEnd: ae, Strand: data[7], BlastzScore: bs, RefSeq: seq1, QuerySeq: seq2}
		answer = append(answer, curr)
	}
	return answer, nil
}

func AxtToSnp(axtFile *Axt) []*vcf.Snp {
	var answer []*vcf.Snp
	var curr *vcf.Snp
	rCount := axtFile.RefStart - 1
	qCount := axtFile.AligningStart - 1
	//rLastMatch := rCount
	//qLastMatch := qCount
	for i := 0; i < len(axtFile.RefSeq); i++ {
		if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), "-") != 0 && strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i])), "-") != 0 {
			rCount++
			qCount++
			//snp mismatch
			if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), strings.ToUpper(string(axtFile.QuerySeq[i]))) != 0 {
				curr = &vcf.Snp{RefName: axtFile.RefName, RefPos: rCount, RefSub: strings.ToUpper(string(axtFile.RefSeq[i])), QueryName: axtFile.AligningName, QueryPos: qCount, QuerySub: strings.ToUpper(string(axtFile.QuerySeq[i]))}
				answer = append(answer, curr)
			}
		}
		//deletion in reference
		if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), "-") == 0 {
			qCount++
			curr = &vcf.Snp{RefName: axtFile.RefName, RefPos: rCount, RefSub: strings.ToUpper(string(axtFile.RefSeq[i])), QueryName: axtFile.AligningName, QueryPos: qCount, QuerySub: strings.ToUpper(string(axtFile.QuerySeq[i]))}
			answer = append(answer, curr)
		}
		//insertion in reference
		if strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i])), "-") == 0 {
			rCount++
			curr = &vcf.Snp{RefName: axtFile.RefName, RefPos: rCount, RefSub: strings.ToUpper(string(axtFile.RefSeq[i])), QueryName: axtFile.AligningName, QueryPos: qCount, QuerySub: strings.ToUpper(string(axtFile.QuerySeq[i]))}
			answer = append(answer, curr)
			
		}
	}
	return answer	
}

func AxtToVcf(axtFile *Axt) []*vcf.Vcf {
	var answer []*vcf.Vcf
	var curr *vcf.Vcf
	rCount := axtFile.RefStart - 1
	qCount := axtFile.AligningStart - 1
	
	//rLastMatch := rCount
	//qLastMatch := qCount
	for i := 0; i < len(axtFile.RefSeq); i++ {
		var infoTag string
		if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), "-") != 0 && strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i])), "-") != 0 {
			rCount++
			qCount++
			//snp mismatch
			if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), strings.ToUpper(string(axtFile.QuerySeq[i]))) != 0 {
				infoTag = "POS=" + strconv.FormatInt(qCount, 10)
				curr = &vcf.Vcf{Chr: axtFile.RefName, Pos: rCount, Id: axtFile.AligningName, Ref: strings.ToUpper(string(axtFile.RefSeq[i])), Alt: strings.ToUpper(string(axtFile.QuerySeq[i])), Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=SNP", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
				//fmt.Println(snps[i].RefSub, snps[i].QuerySub)
				answer = append(answer, curr)
			}
		}
		//insertion in VCF record
		if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i])), "-") == 0 {
			var altTmp string
			qCount++
			//var refTmp string
			altTmp = strings.ToUpper(string(axtFile.RefSeq[i-1]))
			for j := 0; j < len(axtFile.RefSeq); j++ {
				if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i+j])), "-") == 0 && strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i+j])), "-") != 0 {
					
					altTmp = altTmp + strings.ToUpper(string(axtFile.QuerySeq[i+j]))
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RefName, Pos: rCount, Id: axtFile.AligningName, Ref: strings.ToUpper(string(axtFile.RefSeq[i-1])), Alt: altTmp, Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=INS", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
					answer = append(answer, curr)
					i = i + j
					break
				}
			}
		}
		//deleteion vcf record
		if strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i])), "-") == 0 {
			//var refTmp string
			var altTmp string
			rCount++
			altTmp = strings.ToUpper(string(axtFile.RefSeq[i]))
			for j := 0; j < len(axtFile.RefSeq); j++ {
				if strings.Compare(strings.ToUpper(string(axtFile.RefSeq[i+j])), "-") != 0 && strings.Compare(strings.ToUpper(string(axtFile.QuerySeq[i+j])), "-") == 0 {
					rCount++
					altTmp = altTmp + strings.ToUpper(string(axtFile.RefSeq[i+j]))
				} else {
					curr = &vcf.Vcf{Chr: axtFile.RefName, Pos: rCount, Id: axtFile.AligningName, Ref: altTmp, Alt: strings.ToUpper(string(axtFile.RefSeq[i])), Qual: 0, Filter: "PASS", Info: infoTag, Format: "SVTYPE=DEL", Unknown: "GT:DP:AD:RO:QR:AO:QA:GL"}
					answer = append(answer, curr)
					//rCount = rCount + int64(j)
					i = i + j
					break
				}

			//altTmp = strings.ToUpper(string(axtFile.RefSeq[rCount])) + strings.ToUpper(string(axtFile.QuerySeq[qCount]))
			
			
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

func CallSnps(axtList []*Axt) []*vcf.Snp {
	var answer []*vcf.Snp
	var curr []*vcf.Snp
	for i := 0; i < len(axtList); i++ {
		curr = AxtToSnp(axtList[i])
		for _, j := range curr {
			answer = append(answer, j)
		}
	}
	return answer
}

func ReadDictionary(filename string, ref map[string]int) ([]*Axt, error) {
	var answer []*Axt
	var curr *Axt
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	var prefix bool
	rline, prefix, err2 = reader.ReadLine()
	for ; err2 != io.EOF; rline, prefix, err2 = reader.ReadLine() {
		line = string(rline[:])
		//header data:
		for strings.HasPrefix(line, "#") {
			rline, prefix, err2 = reader.ReadLine()
			line = string(rline[:])
		}
		data := strings.Split(line, " ")
		var an int64
		var rs int64
		var re int64
		var as int64
		var ae int64
		var bs int64
		var seq1 string
		var seq2 string

		if len(data) != 9 {
			fmt.Println("something went wrong")
		}
		an, _ = strconv.ParseInt(data[0], 10, 64)
		rs, _ = strconv.ParseInt(data[2], 10, 64)
		re, _ = strconv.ParseInt(data[3], 10, 64)
		as, _ = strconv.ParseInt(data[5], 10, 64)
		ae, _ = strconv.ParseInt(data[6], 10, 64)
		bs, _ = strconv.ParseInt(data[8], 10, 64)

		rline, prefix, err2 = reader.ReadLine()
		seq1 = string(rline[:])
		for prefix != false {
			rline, prefix, err2 = reader.ReadLine()
			if err2 != nil {
				log.Fatal(err2)
			}
			seq1 += string(rline[:])
		}

		rline, prefix, err2 = reader.ReadLine()
		seq2 = string(rline[:])
		for prefix != false {
			rline, prefix, err2 = reader.ReadLine()
			if err2 != nil {
				log.Fatal(err2)
			}
			line = string(rline[:])
			seq2 += line
		}

		if len(seq1) != len(seq2) {
			fmt.Println("Sequences are not the same length")
			log.Fatal(err2)
		}

		rline, prefix, err2 = reader.ReadLine()
		blank := string(rline[:])
		if len(blank) != 0 {
			fmt.Println("not a blank line")
		}

		if strings.Compare(data[7], "-") == 0 {
			refLength := int64(ref[data[1]])
			startIndexConvert := refLength - rs + 1
			endIndexConvert := refLength - re + 1
			curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: startIndexConvert, AligningEnd: endIndexConvert, Strand: "+", BlastzScore: bs, RefSeq: seq1, QuerySeq: seq2}
			answer = append(answer, curr)
		} else {
			curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: as, AligningEnd: ae, Strand: data[7], BlastzScore: bs, RefSeq: seq1, QuerySeq: seq2}
			answer = append(answer, curr)
		}
	}
	return answer, nil
}

func FishMap(ref []*fasta.Fasta) map[string]int {
	m := make(map[string]int)
	//var answer []*fasta.Fasta
	var curr *fasta.Fasta
	for i := 0; i < len(ref); i++ {
		curr = ref[i]
		_, ok := m[curr.Name]
		if !ok {
			m[curr.Name] = len(curr.Seq)
		}
	}
	return m
}

func AutomateMergeAxt(files []string) []*Axt {
	var answer []*Axt
	for i := 0; i < len(files); i++ {
		a, _ := ReadIn(files[i])
		for j := range a {
			answer = append(answer, a[j])
		}
	}
	return answer
}

func FilterDir(dir, suffix string) ([]string, error) {
	files, err := ioutil.ReadDir(dir)
	var answer []string
	if err != nil {
		log.Fatal(err)
	}
	for _, f := range files {
		if !f.IsDir() && strings.HasSuffix(f.Name(), suffix) {
			answer = append(answer, f.Name())

		}
	}
	return answer, err
}

func ManualMergeAxt() []*Axt {
	var answer []*Axt
	var listFiles []string
	listFiles = append(listFiles, "chrI_aligned.axt")
	listFiles = append(listFiles, "chrII_aligned.axt")
	listFiles = append(listFiles, "chrIII_aligned.axt")
	listFiles = append(listFiles, "chrIX_aligned.axt")
	listFiles = append(listFiles, "chrM_aligned.axt")
	listFiles = append(listFiles, "chrV_aligned.axt")
	listFiles = append(listFiles, "chrVI_aligned.axt")
	listFiles = append(listFiles, "chrVII_aligned.axt")
	listFiles = append(listFiles, "chrVIII_aligned.axt")
	listFiles = append(listFiles, "chrX_aligned.axt")
	listFiles = append(listFiles, "chrXI_aligned.axt")
	listFiles = append(listFiles, "chrXII_aligned.axt")
	listFiles = append(listFiles, "chrXIII_aligned.axt")
	listFiles = append(listFiles, "chrXIV_aligned.axt")
	listFiles = append(listFiles, "chrXIX_aligned.axt")
	listFiles = append(listFiles, "chrXV_aligned.axt")
	listFiles = append(listFiles, "chrXVI_aligned.axt")
	listFiles = append(listFiles, "chrXVII_aligned.axt")
	listFiles = append(listFiles, "chrXVIII_aligned.axt")
	listFiles = append(listFiles, "chrXX_aligned.axt")
	listFiles = append(listFiles, "chrXXI_aligned.axt")
	for i := 0; i < len(listFiles); i++ {
		a, _ := ReadIn(listFiles[i])
		for j := range a {
			answer = append(answer, a[j])
		}
	}
	return answer
}

//will eventually move to the compare.go

func CompareCoord(alpha *Axt, beta *Axt) int {
	if alpha.RefStart < beta.RefStart {
		return -1
	}
	if alpha.RefStart > beta.RefStart {
		return 1
	}
	return 0
}

func CompareName(alpha string, beta string) int {
	return strings.Compare(alpha, beta)
}

func CompareAxt(alpha *Axt, beta *Axt) int {
	compareStorage := CompareName(alpha.RefName, beta.RefName)
	if compareStorage != 0 {
		return compareStorage
	} else {
		return CompareCoord(alpha, beta)
	}
}

func Sort(axts []*Axt) {
	sort.Slice(axts, func(i, j int) bool { return CompareAxt(axts[i], axts[j]) == -1 })
}

func SearchBlastz(query []*Axt) []*Axt {
	m := make(map[string]*Axt)
	var answer []*Axt
	var curr *Axt
	for i := 0; i < len(query); i++ {
		curr = query[i]
		_, ok := m[curr.AligningName]
		if !ok {
			m[curr.AligningName] = curr
		} else {
			if curr.BlastzScore > m[curr.AligningName].BlastzScore {
				m[curr.AligningName] = curr
			}
		}
	}
	for k := range m {
		answer = append(answer, m[k])
	}
	return answer
}

func GrepAxtChr(axtFile []*Axt, chr string) []*Axt {
	var answer []*Axt
	var curr *Axt
	for i := 0; i < len(axtFile); i++ {
		if CompareName(axtFile[i].RefName, chr) == 0 {
			curr = axtFile[i]
			answer = append(answer, curr)
		}
	}
	return answer
}

func OrderContigsFromAxt(fastaFile []*fasta.Fasta, axtFile []*Axt) []*fasta.Fasta {
	var answer []*fasta.Fasta
	var curr *fasta.Fasta
	for i := 0; i < len(axtFile); i++ {
		for j := 0; j < len(fastaFile); j++ {
			if CompareName(axtFile[i].AligningName, fastaFile[j].Name) == 0 {
				curr = fastaFile[j]
				answer = append(answer, curr)
			}
		}
	}
	return answer
}

func AddUnalignedContigs(sorted []*fasta.Fasta, draftUn []*fasta.Fasta) []*fasta.Fasta {
	answer := sorted
	var curr *fasta.Fasta
	m := make(map[string]*fasta.Fasta)
	unAlign := make(map[string]*fasta.Fasta)
	//add sorted fastas as keys
	for i := 0; i < len(sorted); i++ {
		curr = sorted[i]
		m[curr.Name] = curr
	}
	for j := 0; j < len(draftUn); j++ {
		curr = draftUn[j]
		_, ok := m[curr.Name]
		if !ok {
			unAlign[curr.Name] = curr
		}
	}
	for k := range unAlign {
		answer = append(answer, unAlign[k])
	}
	return answer
}

//Takes in a sorted fasta (subset) and the entire set and sorts out the unalained
func UnmatchedContigs(sorted []*fasta.Fasta, draftUn []*fasta.Fasta) []*fasta.Fasta {
	var answer []*fasta.Fasta
	var curr *fasta.Fasta
	m := make(map[string]*fasta.Fasta)
	unAlign := make(map[string]*fasta.Fasta)
	//add sorted fastas as keys
	for i := 0; i < len(sorted); i++ {
		curr = sorted[i]
		m[curr.Name] = curr
	}
	for j := 0; j < len(draftUn); j++ {
		curr = draftUn[j]
		_, ok := m[curr.Name]
		if !ok {
			unAlign[curr.Name] = curr
		}
	}
	for k := range unAlign {
		answer = append(answer, unAlign[k])
	}
	return answer
}

func PrintAxt(input []*Axt) {
	for i := range input {
		//fmt.Printf("%v\t%s\t%v\t%v\t%s\t%v\t%v\t%s\t%v\n",input[i].AlignNumber, input[i].RefName, input[i].RefStart, input[i].RefEnd, input[i].AligningName, input[i].AligningStart, input[i].AligningEnd, input[i].Strand, input[i].BlastzScore)
		fmt.Println(input[i].AlignNumber, input[i].RefName, input[i].RefStart, input[i].RefEnd, input[i].AligningName, input[i].AligningStart, input[i].AligningEnd, input[i].Strand, input[i].BlastzScore)
	}
}

