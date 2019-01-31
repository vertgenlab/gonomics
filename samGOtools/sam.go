package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"strings"
	"os"
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
	Alignment  []string
}


func ReadIn(filename string) ([]*string, error){
//func ReadIn(filename string) ([]*Sam, error){
	//var answer []*Sam
	var tmpHeader []*string
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
			fmt.Println(line)
		}
	}
	//return answer, scanner.Err()
	return tmpHeader, scanner.Err()
}
