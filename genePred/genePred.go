package genePred

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type GenePred struct {
	Id         string
	Symbol     string
	Chrom      string
	Strand     bool
	TxStart    int
	TxEnd      int
	CdsStart   int
	CdsEnd     int
	ExonStarts []int
	ExonEnds   []int
	ExonFrames []uint8 //i made this a slice, it should hold the frame for each exon?
	Score      int
}

func Read(filename string) []*GenePred {
	var line string
	var answer []*GenePred
	var doneReading bool = false

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processGenePredLine(line)
		answer = append(answer, current)
	}
	return answer
}

func processGenePredLine(line string) *GenePred {
	current := GenePred{}

	words := strings.Split(line, "\t")
	current.Id = words[0]
	current.Symbol = words[0]
	current.Chrom = words[1]
	if words[2] == "+" {
		current.Strand = true
	} else if words[2] == "-" {
		current.Strand = false
	} else {
		log.Fatal("no strand specified")
	}
	current.TxStart = common.StringToInt(words[3])
	current.TxEnd = common.StringToInt(words[4])
	current.CdsStart = common.StringToInt(words[5])
	current.CdsEnd = common.StringToInt(words[6])
	//exonNumber := len(current.ExonStarts)
	current.ExonStarts = StringToIntSlice(words[7])
	current.ExonEnds = StringToIntSlice(words[8])
	current.ExonFrames = CalcExonFrame(current)
	current.Score = 0

	//if exonNumber != len(current.ExonStarts) {
	//	log.Print(exonNumber)
	//	log.Fatal("exon number does not equal number of start coordinates")
	//}

	if len(current.ExonStarts) != len(current.ExonEnds) {
		log.Fatal("there are not the same number of exon start positions as exon end positions")
	}

	return &current
}

func StringToIntSlice(text string) []int {
	values := strings.Split(text, ",")
	var answer []int = make([]int, len(values))

	for i := 0; i < len(values); i++ {
		answer[i] = common.StringToInt(values[i])
	}
	return answer
}

func CalcExonFrame(gene GenePred) []uint8 {
	exonStarts := gene.ExonStarts
	exonEnds := gene.ExonEnds
	cdsStart := gene.CdsStart
	var length uint8
	var nextExonLength uint8
	var nextExonFrame uint8
	var exonFrames []uint8
	exonFrames = append(exonFrames, 0)

	//for first exon
	length = uint8(exonEnds[0] - cdsStart)
	exonTwoFrame := length % 3
	if exonTwoFrame > 2 {
		log.Fatal("frame is offset by more than 2 positions")
	}
	exonFrames = append(exonFrames, exonTwoFrame)
	//DEBUG: fmt.Print(exonFrames)

	//for all other exons, which depend on the frame being calculated ahead of this step
	for i := 1; i < len(exonEnds); i++ {
		nextExonLength = uint8(exonEnds[i]-exonStarts[i]) + exonFrames[i]
		nextExonFrame = nextExonLength % 3
		if nextExonFrame > 2 {
			log.Fatal("frame is offset by more than 2 positions")
		}
		exonFrames = append(exonFrames, nextExonFrame)
		//DEBUG: fmt.Print(exonFrames)
	}
	return exonFrames
}
