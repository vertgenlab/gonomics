// Command Group: "General Tools"

package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"runtime/trace"
	"sync"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
)

func usage() {
	fmt.Print(
		"intervalOverlap - A tool to find non/overlapping genomic regions\n\n" +
			"Usage:\n" +
			"  intervalOverlap [options] select.file in.file out.file\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func sequelOverlap(options *Settings) chan *queryAnswer {
	selectChan := interval.GoReadToChan(options.SelectFile)

	var intervals []interval.Interval
	for val := range selectChan {
		intervals = append(intervals, val)
	}

	tree := buildTree(intervals, options.Aggregate)

	queryChan := interval.GoReadToChan(options.Input)
	answerChan := make(chan *queryAnswer, 1000)

	var wg sync.WaitGroup
	for i := 0; i < options.Threads; i++ {
		wg.Add(1)
		go queryWorker(tree, queryChan, answerChan, options.Relationship, &wg, options.MergedOutput)
	}

	// Spawn a goroutine that closes answerChan once all queryWorkers have finished
	go func() {
		wg.Wait()
		close(answerChan)
	}()

	return answerChan
}

func main() {
	//var extend *int = flag.Int("extend", 0, "WIP") //TODO
	var nonOverlap *bool = flag.Bool("nonOverlap", false, "Return records that do NOT have any overlap with records in the select file")
	var threads *int = flag.Int("threads", 1, "Number of threads to use.")
	//var percentOverlap *float64 = flag.Float64("percentOverlap", 0, "WIP") //TODO
	//var baseOverlap *int = flag.Int("baseOverlap", 0, "WIP")               //TODO
	var aggregate *bool = flag.Bool("aggregate", false, "Determine overlap based on the sum of overlapping target records rather than individual target records.")
	var relationship *string = flag.String("relationship", "any", "Choose a specific relationships that target and query records must fulfill to be reported. Use --printRelationships for more information.")
	var mergedOutput *bool = flag.Bool("mergedOutput", false, "Print the input line followed by the corresponding select lines in the outfile.")
	//var swapTargetQuery *bool = flag.Bool("swapTargetQuery", false, "WIP") //TODO
	var helpRelationships *bool = flag.Bool("printRelationships", false, "Show a diagram of the valid interval relationships that can be tested for.")
	cpuprof := flag.String("cpuprof", "", "DEBUG: cpu profile output for use with go tool pprof")
	memprof := flag.String("memprof", "", "DEBUG: mem profile output for use with go tool pprof")
	exectrace := flag.String("exectrace", "", "DEBUG: execution trace output for use with go tool trace")
	flag.Parse()

	if *helpRelationships {
		interval.PrintRelationships()
		return
	}

	if *cpuprof != "" {
		cpu, err := os.Create(*cpuprof)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(cpu)
		defer pprof.StopCPUProfile()
	}

	if !interval.TestValidRelationship(*relationship) {
		interval.PrintRelationships()
		log.Fatalln("ERROR: Invalid relationship", *relationship)
	}

	if *mergedOutput && *nonOverlap {
		log.Fatalln("ERROR: Cannot use both mergedOutput and nonOverlap")
	}

	var expectedNumArgs int = 3
	var err error
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
		return
	}

	selectFile, inFile, outFile := flag.Arg(0), flag.Arg(1), flag.Arg(2)
	options := &Settings{
		Input:      inFile,
		Output:     outFile,
		SelectFile: selectFile,
		//Extend:          *extend,
		NonOverlap: *nonOverlap,
		Threads:    *threads,
		//PercentOverlap:  *percentOverlap,
		//BaseOverlap:     *baseOverlap,
		Aggregate:    *aggregate,
		Relationship: *relationship,
		MergedOutput: *mergedOutput,
		//SwapTargetQuery: *swapTargetQuery,
	}

	if *exectrace != "" {
		runtrace := fileio.EasyCreate(*exectrace)
		err = trace.Start(runtrace)
		exception.PanicOnErr(err)
	}

	answerChan := sequelOverlap(options)
	output := fileio.EasyCreate(options.Output)
	writeToFile(answerChan, output, options.MergedOutput, options.NonOverlap)
	err = output.Close()
	exception.PanicOnErr(err)

	if *exectrace != "" {
		trace.Stop()
	}

	if *memprof != "" {
		mem, err := os.Create(*memprof)
		if err != nil {
			log.Fatal(err)
		}
		pprof.WriteHeapProfile(mem)
		mem.Close()
	}
}
