/*package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

type slurmCheckArray struct {
	begin		[] string
	outToCheck	[] string
	end			[] string
}



//input is not a file in this case
func parseTheInput (slurmArrayFancy string) [] slurmCheckArray {
	inFancy := fileio.EasyOpen(slurmArrayFancy)
	var lastzArray [] lastzArrayStruct
	var doneReading bool
	var err			error
	var line 		string

	for line, doneReading = fileio.EasyNextLine(inFancy); !doneReading; line, doneReading = fileio.EasyNextLine(inFancy){
		currentLine := processFancyLastzLine (line)
		lastzArray = append(lastzArray, currentLine)
	}

	err = inFancy.Close()
	exception.PanicOnErr(err)
	return lastzArray
}



func usage () {
	fmt.Print(
		"ZZZ - used to check for completion of SLURM job arrays. Takes in a fancy version of a job array text file. ZZZ")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	slurmArrayFancy := flag.Arg(0)
	parseTheInput(slurmArrayFancy)
}

 */