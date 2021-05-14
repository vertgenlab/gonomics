//Package chain declares the Chain struct that describes an alignment block and contains functions that operates on chain structs.
package chain

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
	"sync"
)

//Chain alignment fields.
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

//BaseStats is a cigar-like info for alignment block: First number is the length/size of bases, then number of target gaps and finally query gaps.
type BaseStats struct {
	Size   int
	TBases int
	QBases int
}

//SeqChain is a data structure that wraps over a chain channel. The struct also includes target and query sequences implemented as a hash, mapping seqeuence names to []dna.Base.
type SeqChain struct {
	Chains chan *Chain
	TSeq   map[string][]dna.Base
	QSeq   map[string][]dna.Base
}

//HeaderComments stores the comment lines at the beginning of chain alignments into a struct.
type HeaderComments struct {
	HashTag []string
}

//Read processes an input file and returns a header stuct and a slice of chain alignments.
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

//ReadToChan processes a chain text file to a channel of chains.
func ReadToChan(file *fileio.EasyReader, data chan<- *Chain, wg *sync.WaitGroup) {
	for curr, done := NextChain(file); !done; curr, done = NextChain(file) {
		data <- curr
	}
	file.Close()
	wg.Done()
}

func GoReadToChan(filename string) (<-chan *Chain, *HeaderComments) {
	file := fileio.EasyOpen(filename)
	var wg sync.WaitGroup
	data := make(chan *Chain)
	header := ReadHeaderComments(file)
	wg.Add(1)
	go ReadToChan(file, data, &wg)

	go func() {
		wg.Wait()
		close(data)
	}()

	return data, header
}

//GoReadSeqChain will wrap a chain file with target and query fasta seqeunces into the SeqChain struct.
func GoReadSeqChain(filename string, target []fasta.Fasta, query []fasta.Fasta) *SeqChain {
	file := fileio.EasyOpen(filename)
	ans := make(chan *Chain)
	var wg sync.WaitGroup
	wg.Add(1)
	go ReadToChan(file, ans, &wg)

	go func() {
		wg.Wait()
		close(ans)
	}()

	return &SeqChain{
		Chains: ans,
		TSeq:   fasta.ToMap(target),
		QSeq:   fasta.ToMap(query),
	}
}

//WriteToFile will process a chain channel and writes the data to a file. Once WriteToFile finishes ranging over the channel, it will call Done() on the waitGroup. WaitGroup must be set up beforehand.
func WriteToFile(filename string, chaining <-chan *Chain, comments *HeaderComments, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if comments != nil {
		WriteHeaderComments(file, comments)
	}
	for data := range chaining {
		WriteChain(file, data)
	}
	wg.Done()
}

func WriteToFileHandle(file io.Writer, rec *Chain) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", ToString(rec))
	common.ExitIfError(err)
}

//Write will write chain slice and any given comments to the top of the file.
func Write(filename string, chaining []*Chain, comments *HeaderComments) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	if comments != nil {
		WriteHeaderComments(file, comments)
	}
	for _, data := range chaining {
		WriteChain(file, data)
	}
}

//WriteHeaderComments will take the HeaderComment struct and write to file.
func WriteHeaderComments(file *fileio.EasyWriter, comments *HeaderComments) {
	var err error
	for _, each := range comments.HashTag {
		_, err = fmt.Fprintf(file, "%s\n", each)
		common.ExitIfError(err)
	}
}

//WriteChain processes a chain struct and writes the data to a file.
func WriteChain(file *fileio.EasyWriter, chaining *Chain) {
	_, err := fmt.Fprintf(file, "%s\n", ToString(chaining))
	common.ExitIfError(err)
}

//ReadHeaderComments will process header comments that sometimes appear at the beginning of chain file and returns a struct.
func ReadHeaderComments(er *fileio.EasyReader) *HeaderComments {
	var line string
	var commments HeaderComments
	for nextBytes, done := er.Peek(1); nextBytes[0] == '#' && done == nil; nextBytes, done = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		commments.HashTag = append(commments.HashTag, line)
	}
	return &commments
}

//ToString will convert a chain struct to original string format.
func ToString(ch *Chain) string {
	var answer string = fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
	//minus one in the loop because last line contains 2 zeros and we do not want to print those
	for i := 0; i < len(ch.Alignment)-1; i++ {
		answer += fmt.Sprintf("%d\t%d\t%d\n", ch.Alignment[i].Size, ch.Alignment[i].TBases, ch.Alignment[i].QBases)
	}
	answer = fmt.Sprintf("%s%d\n", answer, ch.Alignment[len(ch.Alignment)-1].Size)
	return answer
}

//NextChain will read lines in file and return one chain record at a time and a true false determining the EOF.
func NextChain(reader *fileio.EasyReader) (*Chain, bool) {
	header, done := fileio.EasyNextRealLine(reader)
	if done {
		return nil, true
	}
	return NewChain(header, reader), false
}

//NewChain will process text into chain data fields. It will read the first line of the file and assign to header fields and use a reader to read and process the additional lines of the alignment.
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

//chainingHelper is the helper function that will process the chain alignment fields and return the alignment stats.
func chainingHelper(reader *fileio.EasyReader) []*BaseStats {
	var line string
	var data []string
	var answer []*BaseStats
	var curr *BaseStats
	for nextBytes, done := reader.Peek(1); nextBytes[0] != 0 && done == nil; nextBytes, done = reader.Peek(1) {
		line, _ = fileio.EasyNextRealLine(reader)
		data = strings.Split(line, "\t")
		if len(data) == 1 {
			curr = &BaseStats{
				Size:   common.StringToInt(data[0]),
				TBases: 0,
				QBases: 0,
			}
			answer = append(answer, curr)
			//this will advance the reader to the blank line i beliebe the reader will peak at the blank line in the next iteration and exit
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

//printHeader is a pretty print that will return a string containing chain alignment fields without the alignment stats columns.
func printHeader(ch *Chain) string {
	return fmt.Sprintf("chain %d %s %d %c %d %d %s %d %c %d %d %d\n", ch.Score, ch.TName, ch.TSize, common.StrandToRune(ch.TStrand), ch.TStart, ch.TEnd, ch.QName, ch.QSize, common.StrandToRune(ch.QStrand), ch.QStart, ch.QEnd, ch.Id)
}

//Simple swaping of target and query fields
func SwapBoth(ch *Chain) *Chain {
	ch.TName, ch.QName = ch.QName, ch.TName
	ch.TSize, ch.QSize = ch.QSize, ch.TSize
	ch.TStrand, ch.QStrand = ch.QStrand, ch.TStrand
	ch.TStart, ch.QStart = ch.QStart, ch.TStart
	ch.TEnd, ch.QEnd = ch.QEnd, ch.TEnd
	for i := 0; i < len(ch.Alignment); i++ {
		ch.Alignment[i].TBases, ch.Alignment[i].QBases = ch.Alignment[i].QBases, ch.Alignment[i].TBases
	}
	return ch
}

//IsChainFile will check the filename suffix and determine if the file input is a chain.
func IsChainFile(filename string) bool {
	if strings.HasSuffix(filename, ".chain") {
		return true
	} else if strings.HasSuffix(filename, ".chain.gz") {
		return true
	} else {
		return false
	}
}
