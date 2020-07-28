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

// TODO: Move to interval package
func goReadToIntervalChan(inputFile string) chan interval.Interval {
	answer := make(chan interval.Interval)
	go readToIntervalChan(inputFile, answer)
	return answer
}

// TODO: Move to interval package
func readToIntervalChan(inputFile string, send chan interval.Interval) {
	// How the file is read is dependent on the file extension
	filetype := path.Ext(inputFile)

	if filetype == ".gz" {
		// If terminal extension is ".gz" then trim off the gz and get the next extension
		filetype = path.Ext(inputFile[0 : len(inputFile)-len(filetype)])
	}

	switch filetype {
	case ".bed":
		receive := bed.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".axt":
		receive := axt.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".vcf":
		receive, _ := vcf.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}

	case ".sam":
		receive := sam.GoReadToChan(inputFile)
		for val := range receive {
			send <- val
		}
		//TODO: Other filetypes
	}
	close(send)
}
