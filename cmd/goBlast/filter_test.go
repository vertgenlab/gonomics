package main

import (

	"flag"
	"github.com/vertgenlab/gonomics/common"
	"io/ioutil"
	"log"
	"strings"
	"testing"
)

//var readToChanTests []string = []string{"testdata/smallChrI.fa","testdata/smallChrUn.fa"}
var filename string = "/Users/bulbasaur/gonomics/fasta/testdata/gasAcu1.fa"

func TestTranspose(t *testing.T) {
	EfficientTranspose(filename, "testdata/gasAcu1_split", 20)
}

func TestMerge(t *testing.T) {
	flag.Parse()
	var err error
	allFa, err := ioutil.ReadDir("testdata/.")
	common.ExitIfError(err)
	var inList []string
	for _, f := range allFa {
		if strings.HasSuffix(f.Name(), ".fa") {
			log.Printf("%s\n", f.Name())
			inList = append(inList, "testdata/"+f.Name())
		}
	}
	FinalContigsToCanu("None", "testdata/toCanu.fasta", inList)
}
