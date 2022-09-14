package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

func hiCarOverlap(geneList string, hiCarData string, outFile string) {
	var hiCarPairs []string = fileio.Read(hiCarData)
	var genes []string = fileio.Read(geneList)
	var err error
	out := fileio.EasyCreate(outFile)

	hiCarGeneMap := make(map[string]int)

	for i, line := range genes {
		hiCarGeneMap[line] = i
	}
	for _, hiCarCont := range hiCarPairs {
		items := strings.Split(hiCarCont, "\t")
		_, ok := hiCarGeneMap[items[9]]
		_, ok2 := hiCarGeneMap[items[10]]
		if ok || ok2 {
			_, err = fmt.Fprintf(out, hiCarCont+"\n")
			exception.PanicOnErr(err)
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"hiCarOverlap - pull out Hi-Car contacts from a Hi-Car .bedpe file file based on a gene list \n" +
			"Usage:\n" +
			"hiCarOverlap geneList.txt hiCarData.bedpe output.bedpe\n")
}

func main() {
	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	geneList := flag.Arg(0)
	hiCarData := flag.Arg(1)
	outFile := flag.Arg(2)

	hiCarOverlap(geneList, hiCarData, outFile)
}
