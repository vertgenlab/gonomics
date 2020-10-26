package genePred

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

type GenePred struct {
	Id         string
	Symbol     string
	Chrom      string
	Strand     byte
	TxStart    int
	TxEnd      int
	CdsStart   int
	CdsEnd     int
	ExonStarts []int
	ExonEnds   []int
	ExonFrames []int
	Score      int
}

func GenePredToString(g *GenePred) string {
	var answer string

	if g.Strand == '+' {
		answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "+", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g), g.Score)
	} else if g.Strand == '-' {
		answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "-", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g), g.Score)
	} else if g.Strand == '.' {
		answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, ".", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g), g.Score)
	}

	return answer
}

////WriteGenePred writes an input GenePred struct to an os.File with a specified number of GenePred fields.
//func WriteGenePred(file *os.File, input *GenePred, fields int) {
//	var err error
//	_, err = fmt.Fprintf(file, "%s\n", GenePredToString(input, fields))
//	common.ExitIfError(err)
//}

//WriteToFileHandle writes an input GenePred struct with a specified number of fields to an io.Writer
func WriteToFileHandle(file io.Writer, records []*GenePred) {
	for _, rec := range records { //take out if we need writeSliceToFileHandle
		var err error
		_, err = fmt.Fprintf(file, "%s\n", GenePredToString(rec))
		//TODO: fmt.Fprintf is slow
		common.ExitIfError(err)
	}
}

////WriteSliceToFileHandle writes a slice of GenePred structs with a specified number of fields to an io.Writer
//func WriteSliceToFileHandle(file io.Writer, records []*GenePred, fields int) {
//	for _, rec := range records {
//		WriteToFileHandle(file, rec, fields)
//	}
//}

//Write writes a slice of GenePred structs with a specified number of fields to a specified filename.
func Write(filename string, records []*GenePred) {
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteToFileHandle(file, records)
}

func Read(filename string) []*GenePred {
	var line string
	var answer []*GenePred
	var doneReading = false

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
	if words[2] == "+" { // || == "."
		current.Strand = '+'
	} else if words[2] == "-" {
		current.Strand = '-'
	} else if words[2] == "." {
		current.Strand = '.'
	} else {
		log.Fatal("no strand specified")
	}
	current.TxStart = common.StringToInt(words[3])
	current.TxEnd = common.StringToInt(words[4])
	current.CdsStart = common.StringToInt(words[5])
	current.CdsEnd = common.StringToInt(words[6])
	current.ExonStarts = StringToIntSlice(words[7])
	current.ExonEnds = StringToIntSlice(words[8])
	current.ExonFrames = CalcExonFrame(&current)
	current.Score = 0
	exonNumber := len(current.ExonStarts)

	if exonNumber != len(current.ExonStarts) {
		//DEBUG: log.Print(exonNumber)
		log.Fatal("exon number does not equal number of start coordinates")
	}

	if len(current.ExonStarts) != len(current.ExonEnds) {
		log.Fatal("there are not the same number of exon start positions as exon end positions")
	}

	return &current
}

func StringToIntSlice(text string) []int {
	values := strings.Split(text, ",")
	var answer = make([]int, len(values)-1)

	for i := 0; i < len(values)-1; i++ {
		answer[i] = common.StringToInt(values[i])
	}
	return answer
}

func CalcExonFrame(gene *GenePred) []int {
	exonStarts := gene.ExonStarts
	exonEnds := gene.ExonEnds
	cdsStart := gene.CdsStart
	var length int
	var nextExonLength int
	var nextExonFrame int
	var exonFrames []int
	exonFrames = append(exonFrames, 0)
	var answer int

	for i := 0; i < len(exonEnds)-1; i++ { // - 1 compensates for UCSC format with ending comma
		if i == 0 {
			//for first exon
			length = (exonEnds[0] - cdsStart) + 1
			exonTwoFrame := length % 3
			if exonTwoFrame > 2 {
				log.Fatal("frame is offset by more than 2 positions")
			}
			if exonTwoFrame == 0 {
				answer = exonTwoFrame
				exonFrames = append(exonFrames, answer)
			} else {
				answer = 3 - exonTwoFrame
				exonFrames = append(exonFrames, answer)
			}
			//DEBUG: log.Print(exonFrames)
		} else {
			//for all other exons, which depend on the frame being calculated ahead of this step
			nextExonLength = exonEnds[i] - exonStarts[i] + exonFrames[i] + 1
			nextExonFrame = nextExonLength % 3
			if nextExonFrame == 0 {
				answer = nextExonFrame
				exonFrames = append(exonFrames, answer)
			} else {
				answer = 3 - nextExonFrame
				exonFrames = append(exonFrames, answer)
			}
			if nextExonFrame > 2 {
				log.Fatal("frame is offset by more than 2 positions")
			}
			//DEBUG: fmt.Print(exonFrames)
		}
	}
	return exonFrames
}
