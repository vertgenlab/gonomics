package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	//"io"
	"log"
	"strings"
	//"time"
)

type slurmCheckArray struct {
	begin		string // whatever is in front of the squiggles
	outToCheck	string // this is the file to be checked
	checkType	string // this is the code for how we check the output files
	end			string // whatever is after the squiggles
}


//parseTheInput parses each line of the input file
func parseTheInput (slurmArrayFancy string) [] slurmCheckArray {
	inFancy := fileio.EasyOpen(slurmArrayFancy)
	var slurmArray [] slurmCheckArray
	var doneReading bool
	var err			error
	var line 		string

	for line, doneReading = fileio.EasyNextLine(inFancy); !doneReading; line, doneReading = fileio.EasyNextLine(inFancy){

		// fatal if there is an empty line in the file.
		if len(line) == 0 {
			log.Fatal("empty line in file, please remove and rerun.")
		}
		currentLine := processFancySlurmLine(line)
		slurmArray = append(slurmArray, currentLine)
	}

	err = inFancy.Close()
	exception.PanicOnErr(err)
	return slurmArray
}


// processFancySlurmLine is a helper function to parseTheInput
func processFancySlurmLine(lineToProcess string) slurmCheckArray {
	var slurmLineParsed slurmCheckArray
	var parsedFirstSquiggle, parsedSecondSquiggle []string

	//parse on the first squiggle into a slice of strings
	parsedFirstSquiggle = strings.Split(lineToProcess, "{")
	slurmLineParsed.begin = parsedFirstSquiggle[0]

	//parse on the second squiggle into a slice of strings
	parsedSecondSquiggle = strings.Split(parsedFirstSquiggle[1], "}")

	//separate the string into "words" aka fields. White space defined the boundary of each field
	fields := strings.Fields(parsedSecondSquiggle[0])

	//third field will be the type of check to perform
	slurmLineParsed.checkType = fields[2]

	//fourth field is the output file we will be checking
	slurmLineParsed.outToCheck = fields[3]

	slurmLineParsed.end = parsedSecondSquiggle[1]

	return slurmLineParsed
}

//TODO
// next; parse the out in the struct to handle the "line +" ?
// make function for creating the sbatch text file
// make function that communicates with terminal to run the sbatch and array on slurm
// make function that handles the line+ option; how to check that file exists? non empty? end in newline?
//				likely these will all be separate functions.
// make function creates a new fancy array file based on whichever didnt run properly the first time


func usage () {
	fmt.Print(
		"slurmCheck - Used to check for completion of SLURM job arrays. Takes in a 'fancy' version of a job array text file.\n" +
			"Outputs a new fancy version of a job array text file with all jobs that did NOT meet the user specified\n" +
			"check.\n" +
			"The four types of checks:\n" +
			"'exists'\tfile must exist\n" +
			"'exists+'\tfile must exist and be non-empty\n" +
			"'line'\t\tfile must end with a complete line\n" +
			"'line+'\t\tfile must end with a complete line and be non-empty\n" +
			"Usage:\n" +
				"inputFancyJobArrayFile.txt\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 1
	flag.Usage = usage
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	slurmArrayFancy := flag.Arg(0)
	var parsed []slurmCheckArray
	parsed = parseTheInput(slurmArrayFancy)
	begin := parsed[0].begin
	out := parsed[0].outToCheck
	check := parsed[0].checkType
	end := parsed[0].end
	fmt.Printf("begin: %s \n out: %s \n check: %s \n end: %s \n", begin, out, check, end)
}


// Code I'm holding:

/*
	var test string
	asRunes := []rune(parsedSecondSquiggle[0])
	length := len(asRunes)
	test = string(asRunes[10:length])
	fmt.Printf("test is: %s \n", test)
*/

/*
// make function for creating the normal slurm array file

func makeNormalArrayFile (slurmArrayFancy string, slurmArray [] slurmCheckArray) {
	var err error
	arrayName := slurmArrayFancy + "_NormalArray"
	arrayFile := fileio.EasyCreate(arrayName)

	for i := range slurmArray {
		WriteToFileHandle(arrayFile, slurmArray[i])
	}


	err = arrayFile.Close()
	exception.PanicOnErr(err)
}

func WriteToFileHandle(file io.Writer, rec lastzArrayStruct) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}
*/