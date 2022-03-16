package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"runtime"
	"sync"
	"testing"
)

const (
	testSelectFile   string = "testdata/imbalance.bed.gz"
	testInputFile    string = "testdata/imbalance.vcf.gz"
	testRelationship string = "any"
)

func linearOverlap(tree map[string]*interval.IntervalNode) {
	var currVcf vcf.Vcf
	var done bool
	file := fileio.EasyOpen(testInputFile)
	var answer []interval.Interval
	for currVcf, done = vcf.NextVcf(file); !done; currVcf, done = vcf.NextVcf(file) {
		answer = interval.Query(tree, currVcf, testRelationship)
		fmt.Sprint(answer)
	}
}

func concurrentOverlap(tree map[string]*interval.IntervalNode, threads int) {
	queryChan, _ := vcf.GoReadToChan(testInputFile)
	answerChan := make(chan queryAnswer, 1000)

	var wg sync.WaitGroup
	for i := 0; i < threads; i++ {
		wg.Add(1)
		go concurrentQueryWorker(tree, queryChan, answerChan, testRelationship, &wg)
	}

	// Spawn a goroutine that closes answerChan once all queryWorkers have finished
	go func() {
		wg.Wait()
		close(answerChan)
	}()

	for answer := range answerChan {
		fmt.Sprint(answer)
	}
}

func concurrentQueryWorker(tree map[string]*interval.IntervalNode, queryChan <-chan vcf.Vcf, answerChan chan<- queryAnswer, relationship string, wg *sync.WaitGroup) {
	var qa queryAnswer
	for query := range queryChan {
		qa.query = query
		qa.answer = interval.Query(tree, query, relationship)
		answerChan <- qa
	}
	wg.Done()
}

func prepTree() map[string]*interval.IntervalNode {
	selectBeds := bed.Read(testSelectFile)
	intervals := make([]interval.Interval, len(selectBeds))
	for i := range selectBeds {
		intervals[i] = selectBeds[i]
	}
	return interval.BuildTree(intervals)
}

func BenchmarkLinearOverlap(b *testing.B) {
	b.StopTimer()
	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		linearOverlap(tree)
	}
}

func BenchmarkConcurrentOverlapThreads1(b *testing.B) {
	b.StopTimer()
	var threads int = 1
	var cores int = runtime.NumCPU()
	//runtime.GOMAXPROCS(threads)
	if threads > cores {
		log.Printf("WARNING: Requested %d threads but only %d cores are available.", threads, cores)
	}

	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		concurrentOverlap(tree, threads)
	}
}

func BenchmarkConcurrentOverlapThreads2(b *testing.B) {
	b.StopTimer()
	var threads int = 2
	var cores int = runtime.NumCPU()
	//runtime.GOMAXPROCS(threads)
	if threads > cores {
		log.Printf("WARNING: Requested %d threads but only %d cores are available.", threads, cores)
	}

	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		concurrentOverlap(tree, threads)
	}
}

func BenchmarkConcurrentOverlapThreads4(b *testing.B) {
	b.StopTimer()
	var threads int = 4
	var cores int = runtime.NumCPU()
	//runtime.GOMAXPROCS(threads)
	if threads > cores {
		log.Printf("WARNING: Requested %d threads but only %d cores are available.", threads, cores)
	}

	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		concurrentOverlap(tree, threads)
	}
}
func BenchmarkConcurrentOverlapThreads8(b *testing.B) {
	b.StopTimer()
	var threads int = 8
	var cores int = runtime.NumCPU()
	//runtime.GOMAXPROCS(threads)
	if threads > cores {
		log.Printf("WARNING: Requested %d threads but only %d cores are available.", threads, cores)
	}

	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		concurrentOverlap(tree, threads)
	}
}

func BenchmarkConcurrentOverlapThreads16(b *testing.B) {
	b.StopTimer()
	var threads int = 16
	var cores int = runtime.NumCPU()
	runtime.GOMAXPROCS(threads)
	if threads > cores {
		log.Printf("WARNING: Requested %d threads but only %d cores are available.", threads, cores)
	}

	var tree map[string]*interval.IntervalNode = prepTree()
	b.StartTimer()
	for i := 0; i < b.N; i++ {
		concurrentOverlap(tree, threads)
	}
}
