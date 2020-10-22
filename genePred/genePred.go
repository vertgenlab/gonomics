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
	Strand     bool
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
	txStart := &g.TxStart
	txEnd := &g.TxEnd
	cdsStart := &g.CdsStart
	cdsEnd := &g.CdsEnd
	exonStarts := &g.ExonStarts
	exonEnds := &g.ExonEnds
	exonFrames := &g.ExonFrames
	score := &g.Score

	if g.Id != "" && g.Chrom != "" && txStart != nil && txEnd != nil && cdsStart != nil && cdsEnd != nil && exonStarts != nil && exonEnds != nil {
		if g.Strand == true {
			answer = fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Chrom, "+", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g))
		} else {
			answer = fmt.Sprintf("%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Chrom, "-", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g))
		}
	}
	if g.Id != "" && g.Symbol != "" && g.Chrom != "" && txStart != nil && txEnd != nil && cdsStart != nil && cdsEnd != nil && exonStarts != nil && exonEnds != nil {
		if g.Strand == true {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "+", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g))
		} else {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "-", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, CalcExonFrame(g))
		}
	}
	if g.Id != "" && g.Symbol != "" && g.Chrom != "" && txStart != nil && txEnd != nil && cdsStart != nil && cdsEnd != nil && exonStarts != nil && exonEnds != nil && exonFrames != nil {
		if g.Strand == true {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "+", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, g.ExonFrames)
		} else {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "-", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, g.ExonFrames)
		}
	}
	if g.Id != "" && g.Symbol != "" && g.Chrom != "" && txStart != nil && txEnd != nil && cdsStart != nil && cdsEnd != nil && exonStarts != nil && exonEnds != nil && exonFrames != nil && score != nil {
		if g.Strand == true {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "+", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, g.ExonFrames, g.Score)
		} else {
			answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%v\t%v\t%v\n", g.Id, g.Symbol, g.Chrom, "-", g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonStarts, g.ExonEnds, g.ExonFrames, g.Score)
		}
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
	current.ExonStarts = StringToIntSlice(words[7])
	current.ExonEnds = StringToIntSlice(words[8])
	current.ExonFrames = CalcExonFrame(&current)
	current.Score = 0
	exonNumber := len(current.ExonStarts)

	if exonNumber != len(current.ExonStarts) {
		log.Print(exonNumber)
		log.Fatal("exon number does not equal number of start coordinates")
	}

	if len(current.ExonStarts) != len(current.ExonEnds) {
		log.Fatal("there are not the same number of exon start positions as exon end positions")
	}

	return &current
}

func StringToIntSlice(text string) []int {
	values := strings.Split(text, ",")
	var answer = make([]int, len(values)) //add -1 for ucsc genePreds

	for i := 0; i < len(values); i++ {
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

	//for first exon
	length = exonEnds[0] - cdsStart
	exonTwoFrame := length % 3
	if exonTwoFrame > 2 {
		log.Fatal("frame is offset by more than 2 positions")
	}
	exonFrames = append(exonFrames, exonTwoFrame)
	//DEBUG: fmt.Print(exonFrames)

	//for all other exons, which depend on the frame being calculated ahead of this step
	for i := 1; i < len(exonEnds)-1; i++ {
		nextExonLength = exonEnds[i] - exonStarts[i] + exonFrames[i]
		nextExonFrame = nextExonLength % 3
		if nextExonFrame > 2 {
			log.Fatal("frame is offset by more than 2 positions")
		}
		exonFrames = append(exonFrames, nextExonFrame)
		//DEBUG: fmt.Print(exonFrames)
	}
	return exonFrames
}
