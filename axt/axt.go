package axt

import (
	"bufio"
	//"flag"
	"fmt"
	"sort"
	//"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
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
	//RefSeq 		  string
	//QuerySeq      string
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
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, " ")
		//fmt.Println("there is data here")
		//header data:

		var an int64
		var rs int64
		var re int64
		var as int64
		var ae int64
		var bs int64
		if strings.HasPrefix(line, "#") {
			fmt.Println("# header here")

		}
		if len(data) == 9 {
			fmt.Println(data)
			an, _ = strconv.ParseInt(data[0], 10, 64)
			rs, _ = strconv.ParseInt(data[2], 10, 64)
			re, _ = strconv.ParseInt(data[3], 10, 64)
			as, _ = strconv.ParseInt(data[5], 10, 64)
			ae, _ = strconv.ParseInt(data[6], 10, 64)
			bs, _ = strconv.ParseInt(data[8], 10, 64)
			curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: as, AligningEnd: ae, Strand: data[7], BlastzScore: bs}
			answer = append(answer, curr)
		} 
		if strings.HasPrefix(line, "A") || strings.HasPrefix(line, "T") || strings.HasPrefix(line, "C") || strings.HasPrefix(line, "G") {
			fmt.Println(len(line))
		}
		
		

	}
	return answer, nil
}

func SeqTruth(axtLine []string) int {
	if len(axtLine) == 0 && strings.Compare(axtLine[0], "") == 0 {
		return 0
	} else {
		return -1
	}
}
/*
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
	for ; err2 != io.EOF; rline, _, err2 = reader.ReadLine() {
		line = string(rline[:])
		data := strings.Split(line, " ")
		switch {
		case strings.HasPrefix(line, "#"):
			//don't do anything
		case len(data) == 1:
			//these lines are sequences, and we are not recording them
			//fmt.Println(line)
		case len(line) == 0:
			//blank line
			//fmt.Println("found blank")
		case len(data) == 9:
			an, _ := strconv.ParseInt(data[0], 10, 64)
			rs, _ := strconv.ParseInt(data[2], 10, 64)
			re, _ := strconv.ParseInt(data[3], 10, 64)
			as, _ := strconv.ParseInt(data[5], 10, 64)
			ae, _ := strconv.ParseInt(data[6], 10, 64)
			bs, _ := strconv.ParseInt(data[8], 10, 64)

			if strings.Compare(data[7], "-") == 0 {
				refLength := int64(ref[data[1]])
				startIndexConvert := refLength - rs + 1
				endIndexConvert := refLength - re + 1
				curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: startIndexConvert, AligningEnd: endIndexConvert, Strand: data[7], BlastzScore: bs}
				answer = append(answer, curr)
			} else {
				curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: as, AligningEnd: ae, Strand: data[7], BlastzScore: bs}
				answer = append(answer, curr)
			}

		}
	}
	return answer, nil
}
*/


/*
func ReadNew(filename string) ([]*Axt, error) {
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
	var rline1 []byte
	var rline2 []byte
	var rline3 []byte
	var errorLine1 error
	rline1, _, errorLine1 = reader.ReadLine()
	if errorLine1 != io.EOF {
		break
	}
	data := strings.Split(line, " ")
	line = string(rline[:])

	var errorLine2 error
	rline2, _, errorLine2 = reader.ReadLine()
	line = string(rline[:])

	if errorLine2 error != io.EOF {
		break
	}
	var errorLine3 error
	rline3, _, errorLine3 = reader.ReadLine()
	line = string(rline[:])

	if errorLine3 != io.EOF {
		break
	}

	for ;err2 != io.EOF; rline, _, err2 = reader.ReadLine()  {
		line = string(rline[:])
		data := strings.Split(line, " ")
		//fmt.Println("there is data here")

		switch {
		case strings.HasPrefix(line, "#"):
			//don't do anything
			//fmt.Println("found #")
		case strings.Compare(line, "") == 0:
			fmt.Println("Blank found")
		case len(data) == 1:
			fmt.Println(data)
			//these lines are sequences, and we are not recording them
			//fmt.Println("found sequences")
		case SeqTruth(data) == 0:
			fmt.Println("blank")
			//blank line
			//fmt.Println("found blank")
		case len(data) == 9:
			//fmt.Println("found header line")
			an, _ := strconv.ParseInt(data[0], 10, 64)
			rs, _ := strconv.ParseInt(data[2], 10, 64)
			re, _ := strconv.ParseInt(data[3], 10, 64)
			as, _ := strconv.ParseInt(data[5], 10, 64)
			ae, _ := strconv.ParseInt(data[6], 10, 64)
			bs, _ := strconv.ParseInt(data[8], 10, 64)
			curr = &Axt{AlignNumber: an, RefName: data[1], RefStart: rs, RefEnd: re, AligningName: data[4], AligningStart: as, AligningEnd: ae, Strand: data[7], BlastzScore: bs}
			answer = append(answer, curr)
		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer, nil
} */



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
