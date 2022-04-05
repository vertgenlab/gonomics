package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"log"
	//"os/exec"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

type lastzArrayStruct struct {
	lastzPath 			string
	referenceChromPath	string
	queryChromPath		string
	outputFilePath		string
	scoresPath			string
	lastzFormatSpecs	[] string
}
/*
type lastzSpecs struct {
	format	string
	O		int
	E		int
	T		int
	M		int
	K		int
	L		int
	Y		int
}
 */

func processFancyLastzLine(lineToProcess string) lastzArrayStruct {
	var lastzParsed lastzArrayStruct
	parsed := strings.Split(lineToProcess, " ")
	parsed[0] = lastzParsed.lastzPath
	parsed[1] = lastzParsed.referenceChromPath
	parsed[2] = lastzParsed.queryChromPath
	parsed[3] = lastzParsed.outputFilePath
	parsed[4] = lastzParsed.scoresPath

	for i := 5; i < len(parsed) ; i++ {
		lastzParsed.lastzFormatSpecs = append(lastzParsed.lastzFormatSpecs, parsed[i])
	}

	return lastzParsed
}


func parseTheInput (lastzArrayFancy string) [] lastzArrayStruct {
	inFancy := fileio.EasyOpen(lastzArrayFancy)
	lastzArray := make([] lastzArrayStruct, 1)
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

/*
func createTheTextFile()  string {
	pathTextFile := "nil"
	return pathTextFile
}


func createTheOutputForSlurm() string {
	//writeTheSbatchFile
	sbatchFilePath := "nil"
	return sbatchFilePath
}



func runOnSlurm() {

	arg1 := "sbatch"
	arg2 := "pathSbatchScript"

	cmd := exec.Command(arg1, arg2)
	err := cmd.Run()

	if err != nil {
		log.Fatal(err)
}
}

*/
func createSbatchScript(parsed []lastzArrayStruct, lastzArrayFancy string) {
	var err error
	sbatchScriptName := lastzArrayFancy + "SbatchScript"
	file := fileio.EasyCreate(sbatchScriptName)


	for i := range parsed {
		WriteToFileHandle(file, parsed[i])
	}


	err = file.Close()
	exception.PanicOnErr(err)
}

func WriteToFileHandle(file io.Writer, rec lastzArrayStruct) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}


/*

// WriteToFileHandle writes an input Bed struct to an io.Writer
func WriteToFileHandle(file io.Writer, rec Bed) {
	var err error
	_, err = fmt.Fprintf(file, "%s\n", rec)
	exception.PanicOnErr(err)
}

//Write writes a slice of Bed structs to a specified filename.
func Write(filename string, records []Bed) {
	var err error
	file := fileio.EasyCreate(filename)

	for i := range records {
		WriteToFileHandle(file, records[i])
	}
	err = file.Close()
	exception.PanicOnErr(err)
}

*/
func usage () {
	fmt.Print(
		"ZZZ - used to check for completion of SLURM job arrays. Takes in a fancy version of a job array text file. ZZZ")
	flag.PrintDefaults()
}


// expecting a file from Christi at the moment.
func main() {
	var expectedNumArgs = 1
	flag.Usage = usage
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	lastzArrayFancy := flag.Arg(0)
	parsed := parseTheInput(lastzArrayFancy)
	spec := parsed[1].lastzFormatSpecs[0]
	fmt.Print("specs: %s", spec)
	createSbatchScript(parsed, lastzArrayFancy)
}