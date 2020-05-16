package chain

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"sync"
)

type Chain struct {
	Score     int
	TName     string
	TSize     int
	TStrand   bool
	TStart    int
	TEnd      int
	QName     string
	QSize     int
	QStrand   bool
	QStart    int
	QEnd      int
	Alignment []*DiffBases
	Id        int
}

type DiffBases struct {
	Size  int
	TDiff int
	QDiff int
}

type FileComments struct {
	HashTag []string
}

func Read(filename string) ([]*Chain, *FileComments) {
	file := fileio.EasyOpen(filename)
	var answer []*Chain
	defer file.Close()
	notes := SaveComments(file)
	for block, done := NextChain(file); !done; block, done = NextChain(file) {
		answer = append(answer, block)
	}
	return answer, notes
}

func ReadToChan(reader *fileio.EasyReader, answer chan<- *Chain) {
	for data, err := NextChain(reader); !err; data, err = NextChain(reader) {
		answer <- data
	}
	close(answer)
}

func ChanToFile(filename string, chaining <-chan *Chain, comments *FileComments, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if comments != nil {
		WriteComments(file, comments)
	}
	for data := range chaining {
		WriteChain(file, data)
	}
	wg.Done()
}

func Write(filename string, chaining []*Chain, comments *FileComments) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if comments != nil {
		WriteComments(file, comments)
	}
	for _, data := range chaining {
		WriteChain(file, data)
	}
}

func WriteComments(file *fileio.EasyWriter, comments *FileComments) {
	for _, each := range comments.HashTag {
		_, err := fmt.Fprintf(file, "%s\n", each)
		common.ExitIfError(err)
	}
}

func WriteChain(file *fileio.EasyWriter, chaining *Chain) {
	_, err := fmt.Fprintf(file, "%s\n", ChainToString(chaining))
	common.ExitIfError(err)
}

func SaveComments(er *fileio.EasyReader) *FileComments {
	var line string
	var err error
	var nextBytes []byte
	var commments FileComments
	for nextBytes, err = er.Peek(1); nextBytes[0] == '#' && err == nil; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		commments.HashTag = append(commments.HashTag, line)
	}
	return &commments
}

func NextChain(reader *fileio.EasyReader) (*Chain, bool) {
	header, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	answer := NewChain(header)
	answer.Alignment = chainingHelper(reader)
	return answer, false
}

func NewChain(text string) *Chain {
	data := strings.Split(text, " ")
	if len(data) != 13 {
		log.Fatalf("Error: header line needs to contain 13 data fields\n")
	}
	curr := &Chain{
		Score:   common.StringToInt(data[1]),
		TName:   data[2],
		TSize:   common.StringToInt(data[3]),
		TStrand: common.StringToStrand(data[4]),
		TStart:  common.StringToInt(data[5]),
		TEnd:    common.StringToInt(data[6]),
		QName:   data[7],
		QSize:   common.StringToInt(data[8]),
		QStrand: common.StringToStrand(data[9]),
		QStart:  common.StringToInt(data[10]),
		QEnd:    common.StringToInt(data[11]),
		Id:      common.StringToInt(data[12]),
	}
	return curr
}

func chainingHelper(reader *fileio.EasyReader) []*DiffBases {
	var line string
	var data []string
	var answer []*DiffBases
	var curr *DiffBases
	for nextBytes, err := reader.Peek(1); nextBytes[0] != 0 && err == nil; nextBytes, err = reader.Peek(1) {
		line, _ = fileio.EasyNextRealLine(reader)
		data = strings.Split(line, "\t")
		if len(data) == 1 {
			curr = &DiffBases{
				Size:  common.StringToInt(data[0]),
				TDiff: 0,
				QDiff: 0,
			}
			answer = append(answer, curr)
			line, _ = fileio.EasyNextRealLine(reader)
			break
		}
		if len(data) == 3 {
			curr = &DiffBases{
				Size:  common.StringToInt(data[0]),
				TDiff: common.StringToInt(data[1]),
				QDiff: common.StringToInt(data[2]),
			}
			answer = append(answer, curr)
		}
	}
	return answer
}

func ChainToString(ch *Chain) string {
	var answer string = fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
	for i := 0; i < len(ch.Alignment)-1; i++ {
		answer += fmt.Sprintf("%d\t%d\t%d\n", ch.Alignment[i].Size, ch.Alignment[i].TDiff, ch.Alignment[i].QDiff)
	}
	answer = fmt.Sprintf("%s%d\n", answer, ch.Alignment[len(ch.Alignment)-1].Size)
	return answer
}

func printHeader(ch *Chain) string {
	return fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
}
