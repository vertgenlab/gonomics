package main

import (
	"io"
	"sync"

	"github.com/vertgenlab/gonomics/interval"
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

type fileWriter interface {
	WriteToFileHandle(io.Writer)
}

func buildTree(intervals []interval.Interval, aggregate bool) map[string]*interval.IntervalNode {
	if aggregate {
		interval.MergeIntervals(intervals)
	}
	return interval.BuildTree(intervals)
}

func queryWorker(tree map[string]*interval.IntervalNode, queryChan <-chan interval.Interval, answerChan chan<- *queryAnswer, relationship string, wg *sync.WaitGroup, mergedOutput bool) {
	var answer []interval.Interval
	answerTrue := make([]interval.Interval, 1)
	buf := make([]interval.Interval, 1000)
	numSeen := 0
	for query := range queryChan {
		numSeen++
		if mergedOutput {
			answer = interval.Query(tree, query, relationship)
			answerChan <- &queryAnswer{query, answer}
			continue
		}

		// else if no merged output
		if interval.QueryBool(tree, query, relationship, buf) {
			answer = answerTrue
			answerChan <- &queryAnswer{query, answer}
		} else {
			answer = nil
		}
	}
	wg.Done()
}

func writeToFile(answerChan <-chan *queryAnswer, outfile io.Writer, mergedOutput bool, nonoverlap bool) {
	if mergedOutput {
		for val := range answerChan {
			if len(val.answer) != 0 {
				val.query.(fileWriter).WriteToFileHandle(outfile)
				for _, curr := range val.answer {
					curr.(fileWriter).WriteToFileHandle(outfile)
				}
			}
		}
	} else if nonoverlap {
		for val := range answerChan {
			if len(val.answer) == 0 {
				val.query.(fileWriter).WriteToFileHandle(outfile)
			}
		}
	} else {
		for val := range answerChan {
			if len(val.answer) != 0 {
				val.query.(fileWriter).WriteToFileHandle(outfile)
			}
		}
	}
}
