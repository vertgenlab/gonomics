package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"strings"
	"os"
	//will eventually import align package for 
	//"github.com/vertgenlab/gonomics/align"
)

func usage() {
	fmt.Print(
		"Samtools go edition - b\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	samFile := flag.Arg(0)
	//out, _ := ReadIn(samFile)
	ReadIn(samFile)
//	for i := range out {
//		fmt.Println(out[i])
//	}
}
type Sam struct {
	Header []string
	Alignment  []metaData
}

type metaData struct {
	//query template name
	qName string
	bitFlag int64
	refName string
	currPos string
	//mapping quality
	qMap int64
	// will change to align.Cigar
	cigar string
	rNext string
	posNext int64
	tmpLen int64
	seq string
	qualBase string


}


func ReadIn(filename string) ([]*string, error){
//func ReadIn(filename string) ([]*Sam, error){
	//var answer []*Sam
	var tmpHeader []*string
	//where should i put the star?
	//var tmpAlign []string
	//var tmpAlign []*string
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line = scanner.Text()
		if strings.Contains(line, "@") {
			tmpHeader = append(tmpHeader, &line)
			//for debugging
			//fmt.Println(line)
		} else {
			//still need to convert string to int64
			colData := strings.Split(line, "\t")
			metaData{qName: colData[0], bitFrag: fmt.Sprint("%v",colData[1]), refName: colData[2], currPos: colData[3], qMap: fmt.Sprint("%v",colData[4]), cigar: colData[5], rNext: colData[6], posNext: fmt.Sprint("%v",colData[7]), tmpLen: fmt.Sprint("%v",colData[8]), seq: colData[9], qualBase[10]}
			fmt.Println(colData[4])
			//testing to see if reading in is working
			//for i :=0; i < len(tmpAlign); i++ {
			//	fmt.Println(tmpAlign[i], "\n")
			//}
		}
	}
	//return answer, scanner.Err()
	return tmpHeader, scanner.Err()
}
