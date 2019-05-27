package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"
	"io"
	//will eventually import align package for
	//"github.com/vertgenlab/gonomics/align"
)

func usage() {
	fmt.Print(
		"Samtools go edition - b\n")
	flag.PrintDefaults()
}


type Sam struct {
	//Header    []*string
	//Alignment *metaData
	QName   string
	BitFlag int64
	RefName string
	CurrPos string
	//mapping quality
	QMap int64
	// will change to align.Cigar
	Cigar    string
	RNext    string
	PosNext  int64
	TmpLen   int64
	Seq      string
	QualBase string
}

func samToVcf(alignment[]*Sam) {
	var answer []*vcf.Vcf
	for i := range alignment {
		
	}
}


func ReadIn(filename string) ([]*Sam, error) {
var answer []*Sam
	var curr *Sam
	var line string
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(file)
	if err != nil {
		return nil, err
	}
	var err2 error
	var rline []byte
	for ;err2 != io.EOF; rline, _, err2 = reader.ReadLine()  {
		line = string(rline[:])

		if strings.HasPrefix(line, "@") {
			fmt.Println(line)
		}
		data := strings.Split(line, "\t")

		switch {
		case len(data) == 10:
			bf, _ := strconv.ParseInt(data[1], 10, 64)
			qf, _ := strconv.ParseInt(data[4], 10, 64)
			pn, _ := strconv.ParseInt(data[7], 10, 64)
			tl, _ := strconv.ParseInt(data[8], 10, 64)
			//curr = &Sam{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8], Unknown: data[9]}
			curr = &Sam{QName: data[0], BitFlag: bf, RefName: data[2], CurrPos: data[3], QMap: qf, Cigar: data[5], RNext: data[6], PosNext: pn, TmpLen: tl, Seq: data[9], QualBase: " "}
			answer = append(answer, curr)

		case len(data) == 11:
			bf, _ := strconv.ParseInt(data[1], 10, 64)
			qf, _ := strconv.ParseInt(data[4], 10, 64)
			pn, _ := strconv.ParseInt(data[7], 10, 64)
			tl, _ := strconv.ParseInt(data[8], 10, 64)
			//curr = &Sam{Chr: data[0], Pos: position, Id: data[2], Ref: data[3], Alt: data[4], Qual: 0, Filter: data[6], Info: data[7], Format: data[8], Unknown: data[9]}
			curr = &Sam{QName: data[0], BitFlag: bf, RefName: data[2], CurrPos: data[3], QMap: qf, Cigar: data[5], RNext: data[6], PosNext: pn, TmpLen: tl, Seq: data[9], QualBase: data[10]}
			answer = append(answer, curr)
		default:
			//fmt.Println("unexpected line")
		}
	}
	return answer, nil
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
	align, _ := ReadIn(samFile)

	for i := range align {
			fmt.Println(align[i])
		
	}
	//fmt.Println(len(align))
}

//linebyline := metaData{qName: colData[0], bitFlag: bf, refName: colData[2], currPos: colData[3], qMap: qf, cigar: colData[5], rNext: colData[6], posNext: pn, tmpLen: tl, seq: colData[9], qualBase: colData[10]}
//currSam := Sam{Header: tmpHeader, Alignment: &linebyline}