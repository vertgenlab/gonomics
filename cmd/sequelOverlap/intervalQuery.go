package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"path"
	"sync"
)

type Settings struct { // Draft of input setting struct
	Input           string // implemented
	Output          string // implemented
	SelectFile      string // implemented
	Extend          int
	NonOverlap      bool // implemented // not compatible with MergedOutput
	Threads         int  // implemented
	PercentOverlap  float64
	BaseOverlap     int
	Aggregate       bool   // implemented
	Relationship    string // implemented
	MergedOutput    bool   // implemented // not compatible with NonOverlap
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

func queryWorker(tree map[string]*interval.IntervalNode, queryChan <-chan interval.Interval, answerChan chan<- *queryAnswer, relationship string, wg *sync.WaitGroup) {
	var answer []interval.Interval
	for query := range queryChan {
		answer = interval.Query(tree, query, relationship)
		answerChan <- &queryAnswer{query, answer}
	}
	wg.Done()
}

func writeToFile(answerChan <-chan *queryAnswer, outfile *fileio.EasyWriter, mergedOutput bool, nonoverlap bool) {
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
