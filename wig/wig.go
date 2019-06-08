package wig

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
)

type Wig struct {
	StepType string
	Chrom    string
	Start    float64
	Step     float64
	Values   []float64
}

func Read(filename string) ([]Wig, error) {
	var answer []Wig
	var line string
	var currentWig Wig

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		var lineFields []string = strings.Fields(line)

		if lineFields[0] == "variableStep" {
			log.Fatal("Package wig is not compatible with variableStep wigs")
		} else if lineFields[0] == "fixedStep" {
			if (len(answer) != 0){
				answer = append(answer, currentWig) //write last completed wig to outlist
			}

			var lineFields []string = strings.Fields(line)
			if len(lineFields) != 4 {
				log.Fatalf("Invalid number of arguments, expecting 4, received %d\n", len(lineFields))
			}

			currentWig.StepType = "fixedStep"
			var chromList []string = strings.Split(lineFields[1], "=")
			currentWig.Chrom = chromList[1]
			var startList []string = strings.Split(lineFields[2], "=")
			currentWig.Start, err = strconv.ParseFloat(startList[1], 64)
			var stepList []string = strings.Split(lineFields[3], "=")
			currentWig.Step, err = strconv.ParseFloat(stepList[1], 64)
		} else {
			valueFloat, _ := strconv.ParseFloat(line, 64)
			//if err != nil {
			//	return err
			//}
			currentWig.Values = append(currentWig.Values, valueFloat)
		}
	}
	answer = append(answer, currentWig)
	return answer, scanner.Err()
}

func PrintFirst(rec []Wig) {
	if len(rec) == 0 {
		fmt.Println("Empty Wig")
	} else {fmt.Printf("StepType=%s Chrom=%s Start=%v Step=%v\n", rec[0].StepType, rec[0].Chrom,
		rec[0].Start, rec[0].Step)
		for i := range rec[0].Values {
			fmt.Println(rec[0].Values[i])
		}
	}
}

func Write(filename string, rec []Wig) {
	file, err := os.Create(filename)
		if err != nil{
			log.Fatal(err)
		}
	defer file.Close()

	WriteToFileHandle(file, rec)
}

func WriteToFileHandle(file *os.File, rec []Wig) error {
	var err error
	for i := range rec {
		_, err = fmt.Fprintf(file, "%s chrom=%s start=%v step=%v\n", rec[i].StepType, rec[i].Chrom,
			rec[i].Start, rec[i].Step)
		for j := range rec[i].Values {
			_, err = fmt.Fprintf(file, "%v\n", rec[i].Values[j])
		}
	}
	return err
}
