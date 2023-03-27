package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
	"syscall"
)

func makeMap(gnt string) map[string]string {
	mp := make(map[string]string)
	file := fileio.Read(gnt)

	for _, i := range file {
		columns := strings.Split(i, "\t")
		mp[columns[0]] = columns[1]
	}
	return mp
}

func geneIdToName(geneNameTable string, inFile string, outFile string, ncbi bool, ensembl bool) {
	var found bool
	var a, name string
	var err error
	found = false

	in := fileio.Read(inFile)
	out := fileio.EasyCreate(outFile)
	syscall.Chdir("~/go/src/github.com/vertgenlab/gonomics/cmd/geneIdToName/defaults/")
	var mp map[string]string

	if ncbi {
		mp = makeMap("ncbiRefSeqToGeneNameKeyHg38.txt")
	} else if ensembl {
		mp = makeMap("ensIdToNameKey.txt")
	} else {
		mp = makeMap(geneNameTable)
	}

	var z, x int
	for _, i := range in {
		found = false
		columns := strings.Split(i, "\t")
		for j := range columns {
			name, found = mp[columns[j]]
			if found {
				z++
				columns[j] = name
			}
		}
		if !found {
			x++
		}
		a = strings.Join(columns, "\t")
		fileio.WriteToFileHandle(out, a)
	}
	fmt.Printf("geneIDs lifted: %d\n", z)
	fmt.Printf("geneIDs not found: %d\n", x)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("geneIdToName - finds and replaces gene IDs with gene names according to an input table \n" +
		"All files must be tab delimited. \n" +
		"The input table must be in the format geneID -tab- geneName \n" +
		"If a geneID in the input file is not in the inputTable if will remain unchanged in the output\n" +
		"Usage:\n" +
		"geneIdToName inputTable.txt in.txt out.txt\n" +
		"or\n" +
		"geneIdToName -ncbi/ensembl in.txt out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var ncbi *bool = flag.Bool("ncbi", false, "Use the preset NCBI RefSeq ID to Gene Name conversion table (hg38).")
	var ensembl *bool = flag.Bool("ensembl", false, "Use the preset ENSEMBL ID to Gene Name conversion table (hg38).")
	var expectedNumArgs int

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()

	if *ncbi || *ensembl {
		expectedNumArgs = 2
	} else {
		expectedNumArgs = 3
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	}
	var inputTable, in, out string
	if *ncbi || *ensembl {
		in = flag.Arg(0)
		out = flag.Arg(1)
	} else {
		inputTable = flag.Arg(0)
		in = flag.Arg(1)
		out = flag.Arg(2)
	}
	geneIdToName(inputTable, in, out, *ncbi, *ensembl)
}
