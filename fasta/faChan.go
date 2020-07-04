package fasta

import (
	"sync"
	//"io"
	//"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

type FaChannel struct {
	StreamBuf chan *Fasta
	Wg     *sync.WaitGroup
}

func ReadMultiFilesToChan(faChan *FaChannel, files []string) {
	//var curr *fileio.EasyReader
	var fa *Fasta
	var done bool
	for _, each := range files {
		//TODO: Figure out if we should declare os.File/EasyReader outside loop, can we re-use as a pointer?
		curr := fileio.EasyOpen(each)
		defer curr.Close()
		for fa, done = NextFasta(curr); !done; fa, done = NextFasta(curr) {
			faChan.StreamBuf <- fa
		}
	}
	close(faChan.StreamBuf)
}

/*
func (mc *FaPipe) SafeClose() {
	mc.Once.Do(func() {
		close(mc.StdIn)
	})
}*/

func NewFaChannel() *FaChannel {
	var wg sync.WaitGroup
	return &FaChannel{
		StreamBuf: make(chan *Fasta),
		Wg:     &wg,
	}
}

func ReadToChan(filename string, output chan<- *Fasta) {
	file := fileio.EasyOpen(filename)
	defer file.Close()
	for fa, done := NextFasta(file); !done; fa, done = NextFasta(file) {
		output <- fa
	}
	close(output)
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

func WritingChannel(toWrite <-chan *Fasta, wg *sync.WaitGroup, fileOutput string) {
	writer := fileio.EasyCreate(fileOutput)
	defer writer.Close()
	for fa := range toWrite {
		WriteFasta(writer, fa, 50)

	}
	wg.Done()
}

func GoReadToChan(filename string) <-chan *Fasta {
	output := make(chan *Fasta)
	go ReadToChan(filename, output)
	return output
}
