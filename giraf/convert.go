package giraf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strconv"
	"strings"
)

func GirafToString(g *Giraf) string {
	var answer string
	answer += fmt.Sprintf("%s\t%d\t%d\t%c\t%s\t%s\t%d\t%d\t%s\t%s%s", g.QName, g.QStart, g.QEnd, strandToRune(g.PosStrand), PathToString(g.Path), cigar.ToString(g.Aln), g.AlnScore, g.MapQ, dna.BasesToString(g.Seq), Uint8QualToString(g.Qual), NotesToString(g.Notes))
	return answer
}

func stringToGiraf(line string) *Giraf {
	var curr *Giraf
	data := strings.SplitN(line, "\t", 11)
	if len(data) < 10 {
		log.Fatalf("Error: Expecting at least 10 columns, but only found %d on %s", len(data), line)
	} else {
		curr = &Giraf{
			QName:     data[0],
			QStart:    common.StringToInt(data[1]),
			QEnd:      common.StringToInt(data[2]),
			PosStrand: StringToPos(data[3]),
			Path:      FromStringToPath(data[4]),
			Aln:       cigar.FromString(data[5]),
			AlnScore:  common.StringToInt(data[6]),
			MapQ:      uint8(common.StringToInt(data[7])),
			Seq:       dna.StringToBases(data[8]),
			Qual:      ToQualUint8([]rune(data[9])),
			Notes:     FromStringToNotes(data[10])}
	}
	return curr
}

//TODO: Specific txt formats are still up for discussion
func PathToString(p *Path) string {
	nodeString := make([]string, len(p.Nodes))
	for i, v := range p.Nodes {
		nodeString[i] = strconv.Itoa(int(v))
	}
	return fmt.Sprintf("%d:%s:%d", p.TStart, strings.Join(nodeString, ">"), p.TEnd)
}

func FromStringToPath(column string) *Path {
	words := strings.Split(column, ":")
	if len(words) != 3 {
		log.Fatalf("Error: Needs exact 3 values, only found %d", len(words))
	}
	nodes := strings.Split(words[1], ">")
	answer := Path{TStart: common.StringToInt(words[0]), Nodes: make([]uint32, len(nodes)), TEnd: common.StringToInt(words[2])}

	for i, v := range nodes {
		answer.Nodes[i] = common.StringToUint32(v)
	}
	return &answer
}

func NotesToString(notes []Note) string {
	var answer string = ""
	if len(notes) == 0 {
		return answer
	} else {
		for i := 0; i < len(notes); i++ {
			answer += fmt.Sprintf("\t%s", NoteToString(notes[i]))
		}
	}
	return answer
}

//TODO: will move to fastq package
func ToQualUint8(qual []rune) []uint8 {
	var answer []uint8 = make([]uint8, len(qual))
	for i := 0; i < len(qual); i++ {
		answer[i] = uint8(qual[i])
	}
	return answer
}

func ReverseQualUint8Record(qualScore []uint8) {
	for i, j := 0, len(qualScore)-1; i <= j; i, j = i+1, j-1 {
		qualScore[i], qualScore[j] = qualScore[j], qualScore[i]
	}
}

func Uint8QualToString(qual []uint8) string {
	var answer []rune = make([]rune, len(qual))
	for i := 0; i < len(qual); i++ {
		// SAM format uses ascii offset of 33 to make everything start with individual characters
		// without adding 33 you get values like spaces and newlines
		var asciiOffset uint8 = 33
		answer[i] = rune(qual[i] + asciiOffset)
	}
	fmt.Println(answer)
	fmt.Println(string(answer))
	return string(answer)
}

func strandToRune(posStrand bool) rune {
	if posStrand {
		return '+'
	} else {
		return '-'
	}
}

func StringToPos(column string) bool {
	var answer bool
	if strings.Compare(column, "+") == 0 {
		answer = true
	} else if strings.Compare(column, "-") == 0 {
		answer = false
	} else {
		log.Fatalf("Error did not find stand information in this column...\n")
	}
	return answer
}

func NoteToString(txt Note) string {
	var answer string = fmt.Sprintf("%s:%c:%s", txt.Tag, txt.Type, txt.Value)
	return answer
}

func FromStringToNotes(s string) []Note {
	words := strings.Split(s, "\t")
	answer := make([]Note, len(words))
	var text []string
	for i, v := range words {
		text = strings.SplitN(v, ":", 3)
		answer[i] = Note{Tag: text[0], Type: []rune(text[1])[0], Value: text[2]}
	}
	return answer
}

func qualToString(qual []uint8) string {
	answer := make([]string, len(qual))
	for i, v := range qual {
		answer[i] = strconv.Itoa(int(v))
	}
	return strings.Join(answer, ",")
}

func fromStringToQual(s string) []uint8 {
	words := []rune(s)
	answer := make([]uint8, 0)
	for _, v := range words {
		answer = append(answer, uint8(v))
	}
	return answer
}
