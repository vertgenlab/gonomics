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
	outToCheck	string // this will be what's inside the squiggles
	end			string // whatever is after the squiggles
}


//input is a fancy file
func parseTheInput (slurmArrayFancy string) [] slurmCheckArray {
	inFancy := fileio.EasyOpen(slurmArrayFancy)
	var slurmArray [] slurmCheckArray
	var doneReading bool
	var err			error
	var line 		string

	for line, doneReading = fileio.EasyNextLine(inFancy); !doneReading; line, doneReading = fileio.EasyNextLine(inFancy){
		currentLine := processFancySlurmLine (line)
		slurmArray = append(slurmArray, currentLine)
	}

	err = inFancy.Close()
	exception.PanicOnErr(err)
	return slurmArray
}



func processFancySlurmLine(lineToProcess string) slurmCheckArray {
	var slurmLineParsed slurmCheckArray
	//fmt.Printf("lineToProcess is : %s \n", lineToProcess)
	parsedFirstSquiggle := strings.Split(lineToProcess, "{")
	// parsedFirstSquiggle = stuff + out } end
	slurmLineParsed.begin = parsedFirstSquiggle[0]
	parsedSecondSquiggle := strings.Split(parsedFirstSquiggle[1], "}")
	// parsedSecondSquiggle = out + end
	parsedSecondSquiggle[0] = slurmLineParsed.outToCheck
	fmt.Printf("parsedFirstSquiggle[0] is: %s \n", parsedSecondSquiggle[0])
	//parsedSecondSquiggle[1] = slurmLineParsed.end

	return slurmLineParsed
}

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

// make function for creating the sbatch text file
// make function(s) that checks for the line+ in the .out field
// make function that communicates with terminal to run the sbatch and array on slurm
// make function that handles the line+ option; how to check that file exists? non empty? end in newline?
//				likely these will all be separate functions.
// make function that handles the line+ output to create a new array file and new sbatch script file. this can iterate over until the line+
//			output is zero.


func usage () {
	fmt.Print(
		"ZZZ - used to check for completion of SLURM job arrays. Takes in a fancy version of a job array text file. ZZZ")
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
	parsed := parseTheInput(slurmArrayFancy)
	//out := parsed[0].outToCheck
	//fmt.Printf("outFile: %s", out)
	begin := parsed[0].begin
	fmt.Printf("begin: %s", begin)
}

/*
No new lines in the actual lastzWriter output
/hpc/group/vertgenlab/softwareShared/lastz-master/src/lastz
{check in line+ /hpc/group/vertgenlab/vertebrateConservation/pairwise/hg38.byChrom/chr10.fa}
{check in line + /hpc/group/vertgenlab/vertebrateConservation/pairwise/mm39.byChrom/chr5.fa}
--output={check out line+ /hpc/group/vertgenlab/vertebrateConservation/pairwise/hg38.mm39/chr10/chr5.chr10.axt}
--scores=/data/lowelab/RefGenomes/hg38/lastzFiles/human_chimp_v2.mat
--format=axt O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000
 */