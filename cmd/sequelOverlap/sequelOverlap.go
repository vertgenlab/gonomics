package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"sync"
)

func usage() {
	fmt.Print(
		"sequelOverlap - A tool to find non/overlapping genomic regions\n\n" +
			"Usage:\n" +
			"  sequelOverlap [options] select.file in.file out.file\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func sequelOverlap(init *Settings) chan *queryAnswer {
	selectChan := goReadToIntervalChan(init.SelectFile)

	var intervals []interval.Interval
	for val := range selectChan {
		intervals = append(intervals, val)
	}

	tree := buildTree(intervals, init.Aggregate)

	queryChan := goReadToIntervalChan(init.Input)
	answerChan := make(chan *queryAnswer, 1000) // TODO: benchmark buffer size

	var wg sync.WaitGroup
	for i := 0; i < init.Threads; i++ {
		wg.Add(1)
		go queryWorker(tree, queryChan, answerChan, init.Relationship, &wg)
	}

	// Spawn a goroutine that closes answerChan once all queryWorkers have finished
	go func() {
		wg.Wait()
		close(answerChan)
	}()

	return answerChan
}

func main() {
	var extend *int = flag.Int("extend", 0, "WIP") //TODO
	var nonOverlap *bool = flag.Bool("nonOverlap", false, "Return records that do NOT have any overlap with records in the select file")
	var threads *int = flag.Int("threads", 1, "Number of threads to use.")
	var percentOverlap *float64 = flag.Float64("percentOverlap", 0, "WIP") //TODO
	var baseOverlap *int = flag.Int("baseOverlap", 0, "WIP")               //TODO
	var aggregate *bool = flag.Bool("aggregate", false, "Determine overlap based on the sum of overlapping target records rather than individual target records.")
	var relationship *string = flag.String("relationship", "any", "Choose a specific relationships that target and query records must fulfill to be reported. Use --printRelationships for more information.")
	var mergedOutput *bool = flag.Bool("mergedOutput", false, "Print the input line followed by the corresponding select lines in the outfile.")
	var swapTargetQuery *bool = flag.Bool("swapTargetQuery", false, "WIP") //TODO
	var helpRelationships *bool = flag.Bool("printRelationships", false, "Show a diagram of the valid interval relationships that can be tested for.")
	flag.Parse()

	if *helpRelationships {
		interval.PrintRelationships()
		return
	}

	if !interval.TestValidRelationship(*relationship) {
		interval.PrintRelationships()
		log.Fatalln("ERROR: Invalid relationship", *relationship)
	}

	if *mergedOutput && *nonOverlap {
		log.Fatalln("ERROR: Cannot use both mergedOutput and nonOverlap")
	}

	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	} else {
		selectFile, inFile, outFile := flag.Arg(0), flag.Arg(1), flag.Arg(2)
		init := &Settings{
			Input:           inFile,
			Output:          outFile,
			SelectFile:      selectFile,
			Extend:          *extend,
			NonOverlap:      *nonOverlap,
			Threads:         *threads,
			PercentOverlap:  *percentOverlap,
			BaseOverlap:     *baseOverlap,
			Aggregate:       *aggregate,
			Relationship:    *relationship,
			MergedOutput:    *mergedOutput,
			SwapTargetQuery: *swapTargetQuery,
		}
		answerChan := sequelOverlap(init)

		output := fileio.EasyCreate(init.Output)
		writeToFile(answerChan, output, init.MergedOutput, init.NonOverlap)
	}
}

/*
Considering how many options we might have I think we should make this separate from the main function
func overlapSelect(flag string) {
	switch flag {
	case:

	case:
	}
}*/
