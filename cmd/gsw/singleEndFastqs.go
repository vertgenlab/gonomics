package main

import (
	"bytes"
	"io"
	"log"
	"path/filepath"
	"strings"
	"sync"
	"time"

	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
)

func GswToGiraf(ref *genomeGraph.GenomeGraph, readOne string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64) {
	log.Printf("Single end reads detected...\n")
	log.Printf("Indexing the genome...\n\n")
	seedHash := genomeGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup

	fastqPipe := make(chan fastq.FastqBig, 824)
	girafPipe := make(chan giraf.Giraf, 824)
	go readFqGsw(readOne, fastqPipe)
	log.Printf("Scoring matrix used:\n%s\n", genomeGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning %s to genome graph...", strings.Split(filepath.Base(readOne), ".")[0])
	start := time.Now()
	for i := 0; i < threads; i++ {
		go genomeGraph.RoutineFqToGiraf(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, girafPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go writeSingleGiraf(output, girafPipe, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(girafPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func GswToSam(ref *genomeGraph.GenomeGraph, readOne string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64, header sam.Header) {
	log.SetFlags(log.Ldate | log.Ltime)
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := genomeGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	fastqPipe := make(chan fastq.FastqBig, 824)
	samPipe := make(chan sam.Sam, 824)
	go readFqGsw(readOne, fastqPipe)

	log.Printf("Scoring matrix used:\n%s\n", genomeGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning sequence to genome graph...")
	start := time.Now()
	for i := 0; i < threads; i++ {
		go genomeGraph.RoutineGirafToSamSingle(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.WriteFromChan(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func readFqGsw(filename string, answer chan<- fastq.FastqBig) {
	readOne := fileio.NewByteReader(filename)
	for fq, done := fastq.ReadFqBig(readOne); !done; fq, done = fastq.ReadFqBig(readOne) {
		answer <- fq
	}
	close(answer)
}

func writeSingleGiraf(filename string, input <-chan giraf.Giraf, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer
	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	var err error
	for g := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		_, err = buf.WriteString(giraf.ToString(&g))
		common.ExitIfError(err)
		err = buf.WriteByte('\n')
		common.ExitIfError(err)
		io.Copy(file, buf)
		buf.Reset()
		simplePool.Put(buf)
	}
	file.Close()
	wg.Done()
}
