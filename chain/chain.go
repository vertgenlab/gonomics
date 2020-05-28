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
	Alignment []*BaseStats
	Id        int
}

type BaseStats struct {
	Size   int
	TBases int
	QBases int
}

type HeaderComments struct {
	HashTag []string
}

func Read(filename string) ([]*Chain, *HeaderComments) {
	file := fileio.EasyOpen(filename)
	var answer []*Chain
	defer file.Close()
	notes := ReadHeaderComments(file)
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

func WriteToFile(filename string, chaining <-chan *Chain, comments *HeaderComments, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for data := range chaining {
		WriteChain(file, data)
	}
	wg.Done()
}

func Write(filename string, chaining []*Chain, comments *HeaderComments) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for _, data := range chaining {
		WriteChain(file, data)
	}
}

func WriteHeaderComments(file *fileio.EasyWriter, comments *HeaderComments) {
	var err error
	for _, each := range comments.HashTag {
		_, err = fmt.Fprintf(file, "%s\n", each)
		common.ExitIfError(err)
	}
}

func WriteChain(file *fileio.EasyWriter, chaining *Chain) {
	_, err := fmt.Fprintf(file, "%s\n", ToString(chaining))
	common.ExitIfError(err)
}

func ReadHeaderComments(er *fileio.EasyReader) *HeaderComments {
	var line string
	var err error
	var nextBytes []byte
	var commments HeaderComments
	for nextBytes, err = er.Peek(1); nextBytes[0] == '#' && err == nil; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		commments.HashTag = append(commments.HashTag, line)
	}
	return &commments
}

func ToString(ch *Chain) string {
	var answer string = fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
	//minus one in the loop because last line contains 2 zeros and we do not want to print those
	for i := 0; i < len(ch.Alignment)-1; i++ {
		answer += fmt.Sprintf("%d\t%d\t%d\n", ch.Alignment[i].Size, ch.Alignment[i].TBases, ch.Alignment[i].QBases)
	}
	answer = fmt.Sprintf("%s%d\n", answer, ch.Alignment[len(ch.Alignment)-1].Size)
	return answer
}

func NextChain(reader *fileio.EasyReader) (*Chain, bool) {
	header, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	return NewChain(header, reader), false
}

func NewChain(text string, reader *fileio.EasyReader) *Chain {
	data := strings.Split(text, " ")
	if len(data) != 13 {
		log.Fatalf("Error: header line needs to contain 13 data fields\n")
	}
	curr := &Chain{
		Score:     common.StringToInt(data[1]),
		TName:     data[2],
		TSize:     common.StringToInt(data[3]),
		TStrand:   common.StringToStrand(data[4]),
		TStart:    common.StringToInt(data[5]),
		TEnd:      common.StringToInt(data[6]),
		QName:     data[7],
		QSize:     common.StringToInt(data[8]),
		QStrand:   common.StringToStrand(data[9]),
		QStart:    common.StringToInt(data[10]),
		QEnd:      common.StringToInt(data[11]),
		Alignment: chainingHelper(reader),
		Id:        common.StringToInt(data[12]),
	}
	return curr
}

func chainingHelper(reader *fileio.EasyReader) []*BaseStats {
	var line string
	var data []string
	var answer []*BaseStats
	var curr *BaseStats
	for nextBytes, err := reader.Peek(1); nextBytes[0] != 0 && err == nil; nextBytes, err = reader.Peek(1) {
		line, _ = fileio.EasyNextRealLine(reader)
		data = strings.Split(line, "\t")
		if len(data) == 1 {
			curr = &BaseStats{
				Size:   common.StringToInt(data[0]),
				TBases: 0,
				QBases: 0,
			}
			answer = append(answer, curr)
			//this will advance the reader to the blank line
			//i beliebe the reader will peak at the blank line in the next iteration and exit
			line, _ = fileio.EasyNextRealLine(reader)
			return answer
		} else if len(data) == 3 {
			curr = &BaseStats{
				Size:   common.StringToInt(data[0]),
				TBases: common.StringToInt(data[1]),
				QBases: common.StringToInt(data[2]),
			}
			answer = append(answer, curr)
		} else {
			log.Fatalf("Error: expecting alignment data columns to be 3 or 1 but encountered %d\n", len(data))
		}
	}
	return nil
}

func printHeader(ch *Chain) string {
	return fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
}
