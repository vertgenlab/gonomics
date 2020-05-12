package chain

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"sync"
	"strings"
)

type Chain struct {
	Score         int
	TName         string
	TSize         int
	TStrand       bool
	TStart        int
	TEnd          int
	QName         string
	QSize         int
	QStrand       bool
	QStart        int
	QEnd          int
	Data          []*Ungapped
	LastUnGapSize int
	Id            int
}

type Ungapped struct {
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
	//answer := make(chan *Chain)
	for data, err := NextChain(reader); !err; data, err = NextChain(reader) {
		answer <- data
	}
	close(answer)
}

func ChanToFile(filename string, chaining <- chan *Chain, comments *FileComments, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if comments != nil {
		WriteComments(file, comments)
	}
	for data := range chaining {
		WriteChain(file, data)
	}
}

func WriteComments(file *fileio.EasyWriter, comments *FileComments) {
	for _,each := range comments.HashTag {
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

func ChainToString(block *Chain) string {
	var answer string = fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", block.Score, block.TName, block.TSize, common.StrandToRune(block.TStrand), block.TStart, block.TEnd, block.QName, block.QSize, common.StrandToRune(block.QStrand), block.QStart, block.QEnd, block.Id)
	for i := 0; i < len(block.Data); i++ {
		answer += fmt.Sprintf("%d\t%d\t%d\n", block.Data[i].Size, block.Data[i].TDiff, block.Data[i].QDiff)
	}
	answer = fmt.Sprintf("%s%d\n\n", answer, block.LastUnGapSize)
	return answer
}

func NextChain(reader *fileio.EasyReader) (*Chain, bool) {
	header, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	answer := NewChain(header)
	answer.Data, answer.LastUnGapSize = GatherUngapped(reader)
	return answer, false
}

func NewChain(text string) *Chain {
	data := strings.Split(text, " ")
	if len(data) != 13 {
		log.Fatalf("Error: header line needs to contain 13 data fields\n")
	}
	var curr Chain
	curr.Score = common.StringToInt(data[1])
	curr.TName = data[2]
	curr.TSize = common.StringToInt(data[3])
	curr.TStrand = common.StringToStrand(data[4])
	curr.TStart, curr.TEnd = common.StringToInt(data[5]), common.StringToInt(data[6])
	curr.QName = data[7]
	curr.QSize = common.StringToInt(data[8])
	curr.QStrand = common.StringToStrand(data[9])
	curr.QStart, curr.QEnd = common.StringToInt(data[10]), common.StringToInt(data[11])
	curr.Id = common.StringToInt(data[12])
	return &curr
}

func GatherUngapped(reader *fileio.EasyReader) ([]*Ungapped, int) {
	var line string
	var data []string
	var answer []*Ungapped
	var lastBlock int
	for nextBytes, err := reader.Peek(1); nextBytes[0] != 0 && err == nil; nextBytes, err = reader.Peek(1) {
		line, _ = fileio.EasyNextRealLine(reader)
		data = strings.Split(line, "\t")
		curr := Ungapped{}
		if len(data) == 1 {
			lastBlock = common.StringToInt(data[0])
			line, _ = fileio.EasyNextRealLine(reader)
			break
		}
		if len(data) == 3 {
			curr.Size = common.StringToInt(data[0])
			curr.TDiff = common.StringToInt(data[1])
			curr.QDiff = common.StringToInt(data[2])
			answer = append(answer, &curr)
		}
	}
	return answer, lastBlock
}
