// Command Group: "General Tools"

// Produces a binary matrix for accessibility breadth analysis
package main

import (
	"flag"
	"fmt"
	"log"
	"sort"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
)

func intervalSubsetMatrix(unionFile string, fileListFile string, outFile string, fraction bool, markMultipleOverlaps string) {
	var err error
	var i int
	var j interval.Interval
	var currOverlaps []interval.Interval
	var currIntervalChan <-chan interval.Interval
	var unionIntervals = make([]interval.Interval, 0)
	var currLineString string
	var overlapSize int

	unionRecChan := interval.GoReadToChan(unionFile)
	files := fileio.Read(fileListFile)
	out := fileio.EasyCreate(outFile)

	for curr := range unionRecChan {
		unionIntervals = append(unionIntervals, curr)
	}
	unionTree := interval.BuildTree(unionIntervals)

	//this will be our output feature matrix. mat[row] => col, where rows are genomic regions, cols are files.
	var mat = make(map[string][]float64)

	// if in markMultipleOverlaps mode, create output feature matrix for additional multiple overlap output file
	var matMultipleOverlaps = make(map[string][]float64)

	for _, u := range unionIntervals {
		mat[interval.CoordsToString(u)] = make([]float64, len(files)) //initialize all rows to zeros based on the number of files. Map is indexed by coords of regions.

		if markMultipleOverlaps != "" {
			matMultipleOverlaps[interval.CoordsToString(u)] = make([]float64, len(files))
		}

	}

	for i = range files {
		currIntervalChan = interval.GoReadToChan(files[i])
		for j = range currIntervalChan {
			currOverlaps = interval.Query(unionTree, j, "any") // each j (1 bed region from the file of the column) can overlap 0 or 1 or more bed region of the row
			if len(currOverlaps) > 0 {
				for _, k := range currOverlaps { // k >= 1 would mean j overlaps more than 1 bed region of the row
					if fraction {
						overlapSize = interval.OverlapSize(k, j)
						mat[interval.CoordsToString(k)][i] += float64(overlapSize) / float64(interval.IntervalSize(k)) // add number rather than assign number, in case 1 bed region of the row overlaps multiple bed regions in the files of the columns

						// if in markMultipleOverlaps mode, mark whenever there is an overlap hit
						if markMultipleOverlaps != "" && overlapSize > 0 {
							matMultipleOverlaps[interval.CoordsToString(k)][i] += 1
						}

					} else {
						mat[interval.CoordsToString(k)][i] = 1

						// if in markMultipleOverlaps mode, mark whenever there is an overlap hit
						if markMultipleOverlaps != "" && overlapSize > 0 {
							matMultipleOverlaps[interval.CoordsToString(k)][i] += 1
						}
					}
				}
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
	sort.Strings(keys) //just alphanumeric sort, might change to coordinate sort, though order doesn't matter.

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

	// if in markMultipleOverlaps mode, create additional multiple overlap output file
	if markMultipleOverlaps != "" {
		outMultipleOverlaps := fileio.EasyCreate(markMultipleOverlaps)
		_, err = fmt.Fprintln(outMultipleOverlaps, headerString)
		exception.PanicOnErr(err)
		for x := range keys {
			currLineString = keys[x]
			for y := range mat[keys[x]] {
				currLineString = fmt.Sprintf("%s\t%v", currLineString, matMultipleOverlaps[keys[x]][y])
			}
			_, err = fmt.Fprintln(outMultipleOverlaps, currLineString)
			exception.PanicOnErr(err)
		}
		err = outMultipleOverlaps.Close()
		exception.PanicOnErr(err)
	}
}

func usage() {
	fmt.Print("intervalSubsetMatrix - Produces a matrix for accessibility breadth analysis.\n" +
		"Rows correspond to genomic regions. Columns correspond to queried interval files.\n" +
		"Intervals must be unique in union.interval and each file in files.list (each genomic position can only be in at most 1 interval).\n" +
		"Usage:\n" +
		"intervalSubsetMatrix union.interval files.list out.txt\n" +
		"options:\n",
	)
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var fraction *bool = flag.Bool("fraction", false, "fraction mode outputs a fraction matrix (each entry is a fraction, the number of bases of overlap between the bed region of the row and the file of the column divided by the size of the bed region of the row), unlike the default mode, which outputs a binary matrix (each entry is 1 if the bed region of the row overlaps any bed region in the file of the column, and 0 if otherwise")
	var markMultipleOverlaps *string = flag.String("markMultipleOverlaps", "", "markMultipleOverlaps mode creates an additional output file that indicates the number of overlaps each interval in union.interval has. Use this option to provide the multiple overlap output file name")

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

	intervalSubsetMatrix(unionFile, fileListFile, outFile, *fraction, *markMultipleOverlaps)
}
