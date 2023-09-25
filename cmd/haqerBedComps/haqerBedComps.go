package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"strings"
)

type settings struct {
	bedA             string
	bedB             string
	outFile          string
	list             string
	matrixAverage    string
	matrixComponents string
}

type stats struct {
	percA float64
	percB float64
	avg   float64
}

type matrixLine struct {
	name string
	vals []float64
}

func compareTwo(s settings) {
	var out []string
	a := bed.Read(s.bedA)
	b := bed.Read(s.bedB)
	intervalsA := interval.BedSliceToIntervals(a)
	intervalsB := interval.BedSliceToIntervals(b)
	overlapsA, overlapsB, overlapAverage := interval.IntervalSimilarity(intervalsA, intervalsB)
	header := fmt.Sprintf("percent overlaps of %s in %s\tpercent overlaps of %s in %s\tbedSimilarityScore\n", s.bedA, s.bedB, s.bedA, s.bedB)
	data := fmt.Sprintf("%f\t%f\t%f\n", overlapsA, overlapsB, overlapAverage)
	out = append(out, header)
	out = append(out, data)
	fileio.Write(s.outFile, out)
}

func multipleComparisons(s settings) {
	var a, b []bed.Bed
	var intervalsA, intervalsB []interval.Interval
	var aboveDiagonal bool
	var stat stats
	var outMatrix *fileio.EasyWriter
	var allFiles []string = []string{"x"}

	out := fileio.EasyCreate(s.outFile)
	fileio.WriteToFileHandle(out, "A\tB\tpercent overlaps of A in B\tpercent overlaps of B in A\tbedSimilarityScore")

	if s.matrixAverage != "" {
		outMatrix = fileio.EasyCreate(s.matrixAverage)
	}
	if s.matrixComponents != "" {
		outMatrix = fileio.EasyCreate(s.matrixComponents)
	}

	beds := fileio.Read(s.list)
	for _, i := range beds {
		var l matrixLine
		l.name = i
		a = bed.Read(i)
		intervalsA = interval.BedSliceToIntervals(a)
		for _, j := range beds {
			aboveDiagonal = false
			if i == j {
				aboveDiagonal = true
				fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s\t%f\t%f\t%f", i, j, 1.0, 1.0, 1.0))
				if s.matrixComponents != "" || s.matrixAverage != "" {
					l.vals = append(l.vals, 0)
				}
			}
			b = bed.Read(j)
			intervalsB = interval.BedSliceToIntervals(b)
			stat.percA, stat.percB, stat.avg = interval.IntervalSimilarity(intervalsA, intervalsB)
			fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s\t%f\t%f\t%f", i, j, 1.0, 1.0, 1.0))
			switch {
			case s.matrixAverage != "":
				l.vals = append(l.vals, stat.avg)
			case s.matrixComponents != "" && aboveDiagonal:
				l.vals = append(l.vals, stat.percA)
			case s.matrixComponents != "" && !aboveDiagonal:
				l.vals = append(l.vals, stat.percB)
			}
		}
		writeMatrixLine(l, outMatrix)
	}

	for _, i := range beds {
		allFiles = append(allFiles, i)
	}
	tail := strings.Join(allFiles, "\t")
	fileio.WriteToFileHandle(outMatrix, tail)

	err := out.Close()
	exception.PanicOnErr(err)
	err = outMatrix.Close()
	exception.PanicOnErr(err)
}

func writeMatrixLine(line matrixLine, outMatrix *fileio.EasyWriter) {
	
}

func doWork(s settings) {
	if s.list == "" {
		compareTwo(s)
	} else {
		multipleComparisons(s)
	}
}

func usage() {

}

func main() {
	var list *string = flag.String("list", "", "Provide a list of bed files and perform an IntervalSimilarity test on all possible combinations")
	var matrixAverage *string = flag.String("matrix", "", "Provide a filename for a matrix that will have all files from the -list option on the axes. "+
		"The matrix will be populated with IntervalSimilarity metric. Must be used with -list")
	var matrixComponents *string = flag.String("matrix", "", "Provide a filename for a matrix that will have all files from the -list option on the axes. The matrix will be "+
		"populated with the overlap percentage of A in B above the diagonal and B in A below the diagonal. Must be used with -list")

	var expectedNumArgs int

	if *list != "" {
		expectedNumArgs = 1
	} else {
		expectedNumArgs = 3
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	var s settings
	if *list != "" {
		s = settings{
			list:             *list,
			matrixAverage:    *matrixAverage,
			matrixComponents: *matrixComponents,
			outFile:          flag.Arg(0),
		}
	} else {
		s = settings{
			list:             *list,
			matrixAverage:    *matrixAverage,
			matrixComponents: *matrixComponents,
			bedA:             flag.Arg(0),
			bedB:             flag.Arg(1),
			outFile:          flag.Arg(2),
		}
	}
	doWork(s)
}
