package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/ontology"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
	"sort"
	"strconv"
)

func oboHistogram(oboFile string, gafFilesIndex string, outTable string) {
	gafFileNames := fileio.Read(gafFilesIndex)
	oboRecord, _ := obo.Read(oboFile, false)
	var thisGaf []gaf.Gaf
	var ont map[string]*ontology.Ontology
	var genes int
	var id string
	var ids []string
	var answer []int

	ont = ontology.OboToOntology(oboRecord)
	for _, gafFile := range gafFileNames {
		thisGaf, _ = gaf.Read(gafFile)
		ontology.GeneAssignmentsFromGaf(thisGaf, ont)
	}

	for id = range ont {
		if len(ont[id].Genes) != 0 {
			ids = append(ids, id)
		}
	}
	sort.Strings(ids)
	for i := 0; i < len(ids); i++ {
		genes = len(ont[ids[i]].Genes)
		answer = append(answer, genes)
	}

	var answerString []string
	for a := 0; a < len(answer); a++ {
		answerString = append(answerString, strconv.Itoa(answer[a]))
	}
	fileio.Write(outTable, answerString)
}

func usage() {
	fmt.Print(
		"oboHistogram takes an obo file and a file containing a list of gaf files of interest with one file name on " +
			"each line and outputs the number of genes for each ontology in a single column with the intention of making a histogram of the " +
			"distribution of the number of genes assigned to each GO term\n" +
			"Usage:\n" +
			"oboHistogram input.obo gafFiles.txt output.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	oboFile := flag.Arg(0)
	gafFileIndex := flag.Arg(1)
	outTable := flag.Arg(2)

	oboHistogram(oboFile, gafFileIndex, outTable)
}
