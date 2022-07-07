package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"sort"
)

func intervalSubsetMatrix(unionFile string, fileListFile string, outFile string) {
	var err error
	var i int
	var j interval.Interval
	var currOverlaps []interval.Interval
	var currIntervalChan <-chan interval.Interval
	var unionIntervals = make([]interval.Interval, 0)
	var currLineString string

	unionRecChan := interval.GoReadToChan(unionFile)
	files := fileio.Read(fileListFile)
	out := fileio.EasyCreate(outFile)

	for curr := range unionRecChan {
		unionIntervals = append(unionIntervals, curr)
	}
	unionTree := interval.BuildTree(unionIntervals)

	var mat = make(map[string][]int) //this will be our output feature matrix. mat[row] => col, where rows are genomic regions, cols are files.
	for _, u := range unionIntervals {
		mat[interval.CoordsToString(u)] = make([]int, len(files)) //initialize all rows to zeros based on the number of files. Map is indexed by coords of regions.
	}

	for i = range files {
		currIntervalChan = interval.GoReadToChan(files[i])
		for j = range currIntervalChan {
			currOverlaps = interval.Query(unionTree, j, "any")
			if len(currOverlaps) > 0 {
				mat[interval.CoordsToString(currOverlaps[0])][i] = 1
			}
		}
	}
	var headerString = "Region"
	for k := range files {
		headerString = fmt.Sprintf("%s\t%s", headerString, files[k])
	}
	_, err = fmt.Fprintln(out, headerString)
	exception.PanicOnErr(err)

	keys := make([]string, 0, len(mat))
	for k := range mat {
		keys = append(keys, k)
	}
	sort.Strings(keys) //jusst alphanumeric sort, might change to coordinate sort, though order doesn't matter.

	for x := range keys {
		currLineString = keys[x]
		for y := range mat[keys[x]] {
			currLineString = fmt.Sprintf("%s\t%v", currLineString, mat[keys[x]][y])
		}
		_, err = fmt.Fprintln(out, currLineString)
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print("intervalSubsetMatrix - Produces a binary matrix for accessibility breadth analysis.\n" +
		"Rows correspond to genomic regions. Columns correspond to queried interval files.\n" +
		"Usage:\n" +
		"intervalSubsetMatrix union.interval files.list out.txt\n" +
		"options:\n",
	)
	flag.PrintDefaults()
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

	unionFile := flag.Arg(0)
	fileListFile := flag.Arg(1)
	outFile := flag.Arg(2)

	intervalSubsetMatrix(unionFile, fileListFile, outFile)
}
