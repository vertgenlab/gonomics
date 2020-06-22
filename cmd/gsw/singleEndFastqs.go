package main

import (
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/giraf"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"path/filepath"
	"strings"
	"sync"
	"time"
)

func GswToGiraf(ref *simpleGraph.SimpleGraph, readOne string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64) {
	log.Printf("Single end reads detected...\n")
	log.Printf("Indexing the genome...\n\n")
	seedHash := simpleGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup

	fastqPipe := make(chan *fastq.FastqBig, 824)
	girafPipe := make(chan *giraf.Giraf, 824)
	go fastq.ReadBigToChan(readOne, fastqPipe)
	log.Printf("Scoring matrix used:\n%s\n", simpleGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning %s to genome graph...", strings.Split(filepath.Base(readOne), ".")[0])
	start := time.Now()
	for i := 0; i < threads; i++ {
		go simpleGraph.RoutineFqToGiraf(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, girafPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go giraf.GirafChanToFile(output, girafPipe, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(girafPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}

func GswToSam(ref *simpleGraph.SimpleGraph, readOne string, output string, threads int, seedLen int, stepSize int, scoreMatrix [][]int64, header *sam.SamHeader) {
	log.SetFlags(log.Ldate | log.Ltime)
	log.Printf("Paired end reads detected...\n")

	log.Printf("Indexing the genome...\n\n")
	seedHash := simpleGraph.IndexGenomeIntoMap(ref.Nodes, seedLen, stepSize)
	var wgAlign, wgWrite sync.WaitGroup
	fastqPipe := make(chan *fastq.FastqBig, 824)
	samPipe := make(chan *sam.SamAln, 824)
	go fastq.ReadBigToChan(readOne, fastqPipe)

	log.Printf("Scoring matrix used:\n%s\n", simpleGraph.ViewMatrix(scoreMatrix))
	log.Printf("Aligning with the following settings:\n\t\tthreads=%d, seedLen=%d, stepSize=%d\n\n", threads, seedLen, stepSize)
	wgAlign.Add(threads)
	log.Printf("Aligning sequence to genome graph...")
	start := time.Now()
	for i := 0; i < threads; i++ {
		go simpleGraph.RoutineGirafToSamSingle(ref, seedHash, seedLen, stepSize, scoreMatrix, fastqPipe, samPipe, &wgAlign)
	}
	wgWrite.Add(1)
	go sam.SamChanToFile(samPipe, output, header, &wgWrite)
	wgAlign.Wait()
	stop := time.Now()
	close(samPipe)
	wgWrite.Wait()
	log.Printf("GSW aligner finished in %.1f seconds\n", stop.Sub(start).Seconds())
	log.Printf("Enjoy analyzing your data!\n\n--xoxo GG\n")
}
