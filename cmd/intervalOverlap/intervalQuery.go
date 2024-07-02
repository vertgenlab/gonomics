package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"io"
	"sync"
)

type Settings struct { // Draft of input setting struct
	Input      string // implemented
	Output     string // implemented
	SelectFile string // implemented
	//Extend          int
	NonOverlap       bool    // implemented // not compatible with MergedOutput
	Threads          int     // implemented
	ThresholdOverlap float64 //not compatible with NonOverlap or MergedOutput
	//BaseOverlap     int
	Aggregate    bool   // implemented
	Relationship string // implemented
	MergedOutput bool   // implemented // not compatible with NonOverlap
	//SwapTargetQuery bool
}

type queryAnswer struct {
	query  interval.Interval
	answer []interval.Interval
}

type fileWriter interface {
	WriteToFileHandle(io.Writer)
}

func mergeOutput(v queryAnswer) []string {
	var toWrite []string
	var str string
	for i := range v.answer {
		str = fmt.Sprintf("%v\t%v", v.query, v.answer[i])
		toWrite = append(toWrite, str)
	}
	return toWrite
}

func buildTree(intervals []interval.Interval, aggregate bool) map[string]*interval.IntervalNode {
	if aggregate {
		interval.MergeIntervals(intervals)
	}
	return interval.BuildTree(intervals)
}

func queryWorker(tree map[string]*interval.IntervalNode, queryChan <-chan interval.Interval, answerChan chan<- *queryAnswer, relationship string, wg *sync.WaitGroup, mergedOutput bool, thresholdOverlap float64) {
	var answer []interval.Interval
	buf := make([]interval.Interval, 1000)
	numSeen := 0
	for query := range queryChan {
		numSeen++
		answer = getAnswer(query, tree, relationship, mergedOutput, thresholdOverlap, buf)
		answerChan <- &queryAnswer{query, answer}
	}
	wg.Done()
}

func getAnswer(query interval.Interval, tree map[string]*interval.IntervalNode, relationship string, mergedOutput bool, thresholdOverlap float64, buf []interval.Interval) []interval.Interval {
	var answer []interval.Interval

	// special case to run QueryBool for faster processing
	if !mergedOutput && thresholdOverlap == 0 {
		if interval.QueryBool(tree, query, relationship, buf) {
			return make([]interval.Interval, 1)
		} else {
			return nil
		}
	}

	answer = interval.Query(tree, query, relationship)
	answer = passingThreshold(query, answer, thresholdOverlap)

	if answer == nil {
		return nil
	}

	if mergedOutput {
		return answer
	}

	// else no merged output
	return make([]interval.Interval, 1)
}

func passingThreshold(query interval.Interval, answer []interval.Interval, thresholdOverlap float64) []interval.Interval {
	if len(answer) == 0 || thresholdOverlap == 0 {
		return answer
	}
	var ovSize, intSize float64
	for a := range answer {
		ovSize = float64(interval.OverlapSize(answer[a], query))
		intSize = float64(interval.IntervalSize(query))
		if ovSize/intSize >= thresholdOverlap {
			return make([]interval.Interval, 1)
		}
	}
	return nil
}

func writeToFile(answerChan <-chan *queryAnswer, outfile io.Writer, mergedOutput bool, nonOverlap bool, thresholdOverlap float64) {
	var toWrite []string
	if mergedOutput && thresholdOverlap == 0 {
		for val := range answerChan {
			if len(val.answer) != 0 {
				toWrite = mergeOutput(*val)
				for i := range toWrite {
					fileio.WriteToFileHandle(outfile, toWrite[i])
				}
			}
		}
	} else if nonOverlap {
		for val := range answerChan {
			if len(val.answer) == 0 {
				val.query.(fileWriter).WriteToFileHandle(outfile)
			}
		}
	} else if mergedOutput && thresholdOverlap > 0 {
		for val := range answerChan {
			if len(val.answer) != 0 {
				toWrite = mergeOutput(*val)
				for i := range toWrite {
					fileio.WriteToFileHandle(outfile, toWrite[i])
				}
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
