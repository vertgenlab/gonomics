package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/cmd/geneIdToName/presets"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

func makeCustomMap(gnt string) map[string]string {
	mp := make(map[string]string)
	file := fileio.Read(gnt)

	for _, i := range file {
		columns := strings.Split(i, "\t")
		mp[columns[0]] = columns[1]
	}
	return mp
}

func geneIdToName(geneNameTable string, inFile string, outFile string, ncbi bool, ensembl bool, keepMatching bool) {
	var found, write bool
	var a, name string
	var err error

	in := fileio.Read(inFile)
	out := fileio.EasyCreate(outFile)

	var mp map[string]string

	if ncbi {
		mp = presets.MakeNcbiMap()
	} else if ensembl {
		mp = presets.MakeEnsemblMap()
	} else {
		mp = makeCustomMap(geneNameTable)
	}

	var z, x int
	for _, i := range in {
		write, found = false, false
		columns := strings.Split(i, "\t")
		for j := range columns {
			name, found = mp[columns[j]]
			if found {
				z++
				columns[j] = name
				write = true
			}
		}
		if write && keepMatching {
			a = strings.Join(columns, "\t")
			fileio.WriteToFileHandle(out, a)
		} else if !write && keepMatching {
			x++
		} else {
			if !write {
				x++
			}
			a = strings.Join(columns, "\t")
			fileio.WriteToFileHandle(out, a)
		}
	}
	fmt.Printf("geneIDs lifted: %d\n", z)
	fmt.Printf("geneIDs not found: %d\n", x)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("geneIdToName - find and replace gene IDs with gene names according to an input table or preset NCBI RefSeq / ENSEMBL ID conversion tables\n" +
		"All files must be tab delimited.\n" +
		"The input table must be in the format geneID -tab- geneName\n" +
		"By default a geneID in the input file that is not in the inputTable if will remain unchanged in the output\n" +
		"However, non-matching geneIDs can be excluded from the output with the -keepMatching option\n" +
		"If using a pre-set conversion table, an user input table should not be specified\n" +
		"Usage:\n" +
		"geneIdToName [options] inputTable.txt in.txt out.txt\n" +
		"or\n" +
		"geneIdToName [options] -ncbi/ensembl in.txt out.txt\n")
	flag.PrintDefaults()
}

func main() {
	var ncbi *bool = flag.Bool("ncbi", false, "Use the preset NCBI RefSeq ID to Gene Name conversion table (hg38). Curated set only (NM_* or NR_*). Up to date with version GRCh38, March 28, 2023")
	var ensembl *bool = flag.Bool("ensembl", false, "Use the preset ENSEMBL ID to Gene Name conversion table (hg38). (eg: ENSG00....) Up to date with version GRCh38, March 28, 2023")
	var keepMatching *bool = flag.Bool("keepMatching", false, "Only lines that have a correctly lifted gene name will be kept")
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
	geneIdToName(inputTable, in, out, *ncbi, *ensembl, *keepMatching)
}
