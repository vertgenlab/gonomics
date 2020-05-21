package fasta

import (
	"sync"
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

func GoReadToChan(filename string) <-chan *Fasta {
	output := make(chan *Fasta)
	go ReadToChan(filename, output, nil, true)
	return output
}

func ReadToChan(filename string, output chan<- *Fasta, wg *sync.WaitGroup, closeChannel bool) {
	file := fileio.EasyOpen(filename)
	for fa, done := NextFasta(file); !done; fa, done = NextFasta(file) {
		output <- fa
	}
	if closeChannel {
		close(output)
	}
	if wg != nil {
		wg.Done()
	}
}

func NextFasta(reader *fileio.EasyReader) (*Fasta, bool) {
	var fa Fasta
	line, done := fileio.EasyNextLine(reader)
	if done {
		return nil, true
	} else {
		if strings.HasPrefix(line, ">") {
			fa = Fasta{Name: line[1:len(line)], Seq: nextSeq(reader)}
		}
	}
	return &fa, false
}

func nextSeq(reader *fileio.EasyReader) []dna.Base {
	var line string
	var err error
	var nextBytes []byte
	var answer []dna.Base
	for nextBytes, err = reader.Peek(1); len(nextBytes) > 0 && nextBytes[0] != '>' && err == nil; nextBytes, err = reader.Peek(1) {
		line, _ = fileio.EasyNextLine(reader)
		answer = append(answer, dna.StringToBases(line)...)
	}
	return answer
}

func WritingChannel(file *fileio.EasyWriter, output <-chan *Fasta, wg *sync.WaitGroup) {
	for fa := range output {
		WriteHelper(file, fa, 50)
	}
	wg.Done()
}
