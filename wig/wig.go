//package wig provides functions to read, write, and manipulate wig files.
//more information on the WIG file format can be found at https://genome.ucsc.edu/goldenPath/help/wiggle.html
package wig

import (
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"io"
	"log"
	"strings"
)

//Wig stores information on the chromosome location and step properties of Wig data. Individual wig values are stored in the underlying WigValue struct.
type Wig struct {
	StepType string
	Chrom    string
	Start    int64
	Step     int64
	Values   []*WigValue //{position, score}
}

//WigValue is an internal struct of the Wig data structure that contains information on the Wig value at an individual genomic postion.
type WigValue struct {
	Position int64
	Value    float64
}

//Read generates a Wig data structure from an input filename, provided as a string for a WIG format file.
func Read(filename string) []*Wig {
	var answer []*Wig
	var line string
	var currentWig *Wig
	var doneReading bool = false

	file := fileio.EasyOpen(filename)

	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		var lineFields []string = strings.Fields(line)

		if strings.HasPrefix(line, "#") {
			//do nothing, comment line in wig
		} else if lineFields[0] == "variableStep" {
			currentWig = new(Wig)
			answer = append(answer, currentWig)
			currentWig.StepType = "variableStep"

			var lineFields []string = strings.Fields(line)
			var chromList []string = strings.Split(lineFields[1], "=")
			currentWig.Chrom = chromList[1]

		} else if lineFields[0] == "fixedStep" {
			if len(lineFields) != 4 {
				log.Fatalf("Invalid number of arguments, expecting 4, received %d\n", len(lineFields))
			}

			currentWig = new(Wig)
			answer = append(answer, currentWig)

			var lineFields []string = strings.Fields(line)
			var chromList []string = strings.Split(lineFields[1], "=")
			currentWig.Chrom = chromList[1]
			var startList []string = strings.Split(lineFields[2], "=")
			currentWig.Start = common.StringToInt64(startList[1])
			var stepList []string = strings.Split(lineFields[3], "=")
			currentWig.Step = common.StringToInt64(stepList[1])
			currentWig.StepType = "fixedStep"

		} else {
			//check length = 2 inputs is variable, 1 is fixed
			var lineFields []string = strings.Fields(line)
			currentValue := new(WigValue)

			if len(lineFields) == 1 {
				//fixedStep
				currentValue.Value = common.StringToFloat64(lineFields[0])
				currentValue.Position = currentWig.Start + int64(len(currentWig.Values))*currentWig.Step
			} else if len(lineFields) == 2 {
				//variableStep
				currentValue.Position = common.StringToInt64(lineFields[0])
				currentValue.Value = common.StringToFloat64(lineFields[1])
			}

			currentWig.Values = append(currentWig.Values, currentValue)
		}
	}
	return answer
}

//Prints the first record in a Wig struct. Mainly used for debugging.
func PrintFirst(rec []*Wig) {
	if len(rec) == 0 {
		fmt.Println("Empty Wig")
	} else {
		fmt.Printf("StepType=%s Chrom=%s Start=%v Step=%v\n", rec[0].StepType, rec[0].Chrom,
			rec[0].Start, rec[0].Step)
		if rec[0].StepType == "fixedStep" {
			for i := range rec[0].Values {
				fmt.Println(rec[0].Values[i].Value)
			}
		} else if rec[0].StepType == "variableStep" {
			for i := range rec[0].Values {
				fmt.Printf("%d, %f", rec[0].Values[i].Position, rec[0].Values[i].Value)
			}
		}
	}
}

//Write writes a Wig data structure to a WIG format file at the input filename
func Write(filename string, rec []*Wig) {
	file := fileio.EasyCreate(filename)
	defer file.Close()
	for i := range rec {
		log.Printf("Printing wig object: %d\n", i)
		WriteToFileHandle(file, rec[i])
	}
}

//WriteToFileHandle is an helper function for Write that writes the Wig data structure to a file.
func WriteToFileHandle(file io.Writer, rec *Wig) {
	var err error
	if rec.StepType == "fixedStep" {
		fmt.Println("Printing header")
		_, err = fmt.Fprintf(file, "%s chrom=%s start=%d step=%d\n", rec.StepType, rec.Chrom,
			rec.Start, rec.Step)
		common.ExitIfError(err)
	} else if rec.StepType == "variableStep" {
		_, err = fmt.Fprintf(file, "%s chrom=%s", rec.StepType, rec.Chrom)
		common.ExitIfError(err)
	} else {
		log.Fatalf("Invalid step type for wig.")
	}

	for j := range rec.Values {

		if rec.StepType == "fixedStep" {
			_, err = fmt.Fprintf(file, "%f\n", rec.Values[j].Value)
			common.ExitIfError(err)
		} else if rec.StepType == "variableStep" {
			_, err = fmt.Fprintf(file, "%d\t%f\n", rec.Values[j].Position, rec.Values[j].Value)
			common.ExitIfError(err)
		}

	}
}
