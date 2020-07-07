package main

import (
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"sync"
)

type Settings struct { // Draft of input setting struct
	Extend          int
	Format          string
	NonOverlap      bool // implemented // not compatible with MergedOutput
	Threads         int
	Output          string
	PercentOverlap  float64
	BaseOverlap     int
	Aggregate       bool // implemented
	Relationship    string
	MergedOutput    bool // implemented // not compatible with NonOverlap
	SwapTargetQuery bool
}

type queryAnswer struct {
	query  interval.Interval
	answer []interval.Interval
}

func buildTree(intervals []interval.Interval, aggregate bool) map[string]*interval.IntervalNode {
	if aggregate {
		interval.MergeIntervals(intervals)
	}
	return interval.BuildTree(intervals)
}

func queryWorker(tree map[string]*interval.IntervalNode, queryChan <-chan interval.Interval, answerChan chan<- queryAnswer, relationship string, wg *sync.WaitGroup) {
	var answer []interval.Interval
	for query := range queryChan {
		answer = interval.Query(tree, query, relationship)
		answerChan <- queryAnswer{query, answer}
	}
	wg.Done()
}

func writeToFile(answerChan <-chan queryAnswer, outfile *fileio.EasyWriter, mergedOutput bool, nonoverlap bool) {
	if mergedOutput {
		for val := range answerChan {
			if len(val.answer) != 0 {
				val.query.WriteToFileHandle(outfile)
				for _, curr := range val.answer {
					curr.WriteToFileHandle(outfile)
				}
			}
		}
	} else if nonoverlap {
		for val := range answerChan {
			if len(val.answer) == 0 {
				val.query.WriteToFileHandle(outfile)
			}
		}
	} else {
		for val := range answerChan {
			if len(val.answer) != 0 {
				val.query.WriteToFileHandle(outfile)
			}
		}
	}
}
