// genePred is a package for reading, writing, and manipulating genePred structs
package genePred

import (
	"fmt"
	"io"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
)

// GenePred is a data struct which defines information about a gene. It's based on the format listed by UCSC Genome Browser.
type GenePred struct {
	Id         string
	Symbol     string
	Chrom      string
	Strand     byte
	TxStart    int
	TxEnd      int
	CdsStart   int
	CdsEnd     int
	ExonNum    int
	ExonStarts []int
	ExonEnds   []int
	ExonFrames []int
	Score      int
}

// GenePredToString formats the data held in the genePred struct into a string so it can be written to a file
func GenePredToString(g GenePred) string {
	var answer string

	answer = fmt.Sprintf("%s\t%s\t%s\t%s\t%v\t%v\t%v\t%v\t%v\t%s\t%s\t%s\t%v", g.Id, g.Symbol, g.Chrom, string(g.Strand), g.TxStart, g.TxEnd, g.CdsStart, g.CdsEnd, g.ExonNum, SliceIntToString(g.ExonStarts), SliceIntToString(g.ExonEnds), SliceIntToString(CalcExonFrame(g)), g.Score)
	return answer
}

// WriteToFileHandle writes an input GenePred struct with a specified number of fields to an io.Writer
func WriteToFileHandle(file io.Writer, records []GenePred) {
	for _, rec := range records { //take out if we need writeSliceToFileHandle
		var err error
		_, err = fmt.Fprintf(file, "%s\n", GenePredToString(rec))
		//TODO: fmt.Fprintf is slow
		common.ExitIfError(err)
	}
}

// Write writes a slice of GenePred structs with a specified number of fields to a specified filename.
func Write(filename string, records []GenePred) {
	file := fileio.EasyCreate(filename)
	defer file.Close()

	WriteToFileHandle(file, records)
}

// Read reads a file and processes it into a genePred
func Read(filename string) []GenePred {
	var line string
	var answer []GenePred
	var doneReading = false

	file := fileio.EasyOpen(filename)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		current := processGenePredLine(line)
		answer = append(answer, current)
	}
	return answer
}

// processGenePredLine is a helper function that takes the content of a .gp file and processes the strings into values in the genePred struct
func processGenePredLine(line string) GenePred {
	current := GenePred{}

	words := strings.Split(line, "\t")
	current.Id = words[0]
	current.Symbol = words[0]
	current.Chrom = words[1]
	if words[2] == "+" {
		current.Strand = '+'
	} else if words[2] == "." {
		current.Strand = '.'
	} else if words[2] == "-" {
		current.Strand = '-'
	} else {
		log.Fatal("no strand specified")
	}
	current.TxStart = common.StringToInt(words[3])
	current.TxEnd = common.StringToInt(words[4])
	current.CdsStart = common.StringToInt(words[5])
	current.CdsEnd = common.StringToInt(words[6])
	current.ExonNum = common.StringToInt(words[7])
	if !strings.HasSuffix(words[8], ",") {
		log.Fatal("Exon Starts slice doesn't end in empty string.")
	}
	current.ExonStarts = StringToIntSlice(words[8])
	if !strings.HasSuffix(words[9], ",") {
		log.Fatal("Exon Ends slice doesn't end in empty string.")
	}
	current.ExonEnds = StringToIntSlice(words[9])
	current.ExonFrames = CalcExonFrame(current)
	current.Score = 0

	if current.ExonNum != len(current.ExonStarts) {
		//DEBUG: log.Print(exonNumber)
		log.Fatal("exon number does not equal number of start coordinates")
	}

	if len(current.ExonStarts) != len(current.ExonEnds) {
		log.Fatal("there are not the same number of exon start positions as exon end positions")
	}

	return current
}

// StringToIntSlice is used to processes the string that is provided for multiple exon stops, starts and frames into slices of ints for the corresponding genePred fields
func StringToIntSlice(text string) []int {
	values := strings.Split(text, ",")
	var answer = make([]int, len(values)-1)

	for i := 0; i < len(values)-1; i++ {
		answer[i] = common.StringToInt(values[i])
	}
	return answer
}

// CaclExonFame calculates the frame of each exon in a genePred based on the information in the fields for exonStarts, exonStops and cdsStarts
// TODO: May not work with - strand transcripts
func CalcExonFrame(gene GenePred) []int {
	exonStarts := gene.ExonStarts
	exonEnds := gene.ExonEnds
	cdsStart := gene.CdsStart
	var length int
	var nextExonLength int
	var nextExonFrame int
	var exonFrames []int
	exonFrames = append(exonFrames, 0)
	var answer int

	for i := 0; i < len(exonEnds)-1; i++ { // - 1 compensates for not needing to calculate the frame of the exon after the last
		if i == 0 {
			//for first exon
			length = exonEnds[0] - cdsStart
			exonTwoFrame := length % 3
			if exonTwoFrame > 2 {
				log.Fatal("frame is offset by more than 2 positions")
			}
			if exonTwoFrame == 0 {
				answer = exonTwoFrame
			} else {
				answer = 3 - exonTwoFrame
			}
			exonFrames = append(exonFrames, answer)
			//DEBUG: log.Print(exonFrames)
		} else {
			//for all other exons, which depend on the frame being calculated ahead of this step
			nextExonLength = exonEnds[i] - exonStarts[i] - exonFrames[i]
			nextExonFrame = nextExonLength % 3
			if nextExonFrame == 0 {
				answer = nextExonFrame
			} else {
				answer = 3 - nextExonFrame
			}
			exonFrames = append(exonFrames, answer)
			if nextExonFrame > 2 {
				log.Fatal("frame is offset by more than 2 positions")
			}
			//DEBUG: fmt.Print(exonFrames)
		}
	}
	return exonFrames
}

// SliceIntToString is the reverse of StringToIntSlice and is used for writing a slice of ints to a string for writing genePred files
func SliceIntToString(slice []int) string {
	var buffer strings.Builder
	for i := 0; i < len(slice); i++ {
		buffer.WriteString(fmt.Sprintf("%d,", slice[i]))
	}
	return buffer.String()
}
