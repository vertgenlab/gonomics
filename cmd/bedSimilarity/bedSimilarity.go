package main

//TODO : add check for self-overlap bed files
//TODO: add enrichment calculations

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

type matrixLine struct {
	name string
	vals []float64
}

func compareTwo(s settings) {
	var out []string
	a := bed.Read(s.bedA)
	b := bed.Read(s.bedB)
	aName := separatePath(s.bedA)
	bName := separatePath(s.bedB)
	intervalsA := interval.BedSliceToIntervals(a)
	intervalsB := interval.BedSliceToIntervals(b)
	overlapsA, overlapsB, overlapAverage := interval.IntervalSimilarity(intervalsA, intervalsB)
	header := fmt.Sprintf("proportion overlaps of %s in %s\tproportion overlaps of %s in %s\tbedSimilarityScore", aName, bName, bName, aName)
	data := fmt.Sprintf("%f\t%f\t%f", overlapsA, overlapsB, overlapAverage)
	out = append(out, header)
	out = append(out, data)
	fileio.Write(s.outFile, out)
}

func multipleComparisons(s settings) {
	var intervalsA, intervalsB []interval.Interval
	var percA, percB, avg float64
	var outMatrix *fileio.EasyWriter
	var err error
	var j int
	var l matrixLine
	var aName, bName string
	var allFiles []string = []string{"x"}

	out := fileio.EasyCreate(s.outFile)
	fileio.WriteToFileHandle(out, "A\tB\tproportion overlaps of A in B\tproportion overlaps of B in A\tbedSimilarityScore")

	if s.matrixAverage != "" {
		outMatrix = fileio.EasyCreate(s.matrixAverage)
	}
	if s.matrixComponents != "" {
		outMatrix = fileio.EasyCreate(s.matrixComponents)
	}

	files := fileio.Read(s.list)

	if s.matrixAverage != "" || s.matrixComponents != "" {
		for _, i := range files {
			allFiles = append(allFiles, separatePath(i))
		}
		header := strings.Join(allFiles, "\t")
		fileio.WriteToFileHandle(outMatrix, header)
	}

	// Read data from files once and store
	numFiles := len(files)
	intervalData := make([][]interval.Interval, numFiles) // make slice of slices to put interval data in
	names := make([]string, numFiles)                     // make slice of string for file names

	for i, file := range files { // loop through files
		intervalData[i] = interval.BedSliceToIntervals(bed.Read(file)) //put the []interval.Intervals from the files into the correct index in intervalData
		names[i] = separatePath(file)                                  //put the name of the file in the same index as used for the interval data
	}

	for i := range files {
		l.name = names[i]
		aName = names[i]
		l.vals = []float64{}
		intervalsA = intervalData[i]
		for j = range files {
			if files[i] == files[j] {
				if s.matrixComponents != "" || s.matrixAverage != "" {
					l.vals = append(l.vals, 1.0)
				}
				continue
			}
			intervalsB = intervalData[j]
			percA, percB, avg = interval.IntervalSimilarity(intervalsA, intervalsB)
			if j > i {
				bName = names[j]
				fileio.WriteToFileHandle(out, fmt.Sprintf("%s\t%s\t%f\t%f\t%f", aName, bName, percA, percB, avg))
			}

			if s.matrixAverage != "" {
				l.vals = append(l.vals, avg)
			} else if s.matrixComponents != "" {
				l.vals = append(l.vals, percA)
			}
		}
		if s.matrixAverage != "" || s.matrixComponents != "" {
			writeMatrixLine(l, outMatrix)
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
	if s.matrixAverage != "" || s.matrixComponents != "" {
		err = outMatrix.Close()
		exception.PanicOnErr(err)
	}
}

func writeMatrixLine(line matrixLine, outMatrix *fileio.EasyWriter) {
	toWrite := []string{line.name}
	for _, i := range line.vals {
		toWrite = append(toWrite, fmt.Sprintf("%f", i))
	}
	write := strings.Join(toWrite, "\t")
	fileio.WriteToFileHandle(outMatrix, write)
}

func bedSimilarity(s settings) {
	if s.list == "" {
		compareTwo(s)
	} else {
		multipleComparisons(s)
	}
}

func separatePath(path string) string {
	pathSlice := strings.Split(path, "/")
	return pathSlice[len(pathSlice)-1]
}

func usage() {
	fmt.Print("bedSimilarity -- Takes in 2 bed files or a list of bed files and gives similarity statistics based on number of overlaps" +
		"for the input bed files.\n" +
		"If the -list option is used, all combinations of provided files (excluding the comparison with itself) will be performed. " +
		"The statistics are: the proportion of elements that have an overlap an element in the second bed file, " +
		"the proportion of elements in the second bed file that overlap an element in the first bed file, and " +
		"the average of the two overlap percentages (a metric of overall similarity between the two bed files)." +
		"When using the -list option, either matrix option can be used to create a heatmap of similarities scores for all bed files in the list\n\n" +
		"Usage: \n" +
		"haqerBedComps inBedA inBedB outFile\n" +
		"or\n" +
		"haqerBedComps -list list.txt [options] outFile\n\n")
	flag.PrintDefaults()
}

func main() {
	var list *string = flag.String("list", "", "Provide a list of bed files and perform an bedSimilarity test on all possible combinations")
	var matrixAverage *string = flag.String("matrixAverage", "", "Provide a filename for a matrix that will have all files from the -list option on both axes. "+
		"The matrix will be populated with average bedSimilarity metric. Must be used with -list")
	var matrixComponents *string = flag.String("matrixComponents", "", "Provide a filename for a matrix that will have all files from the -list option on both axes. "+
		"The matrix will be populated with the overlap proportion of elements from the file on the vertical axis with the file on the horizontal axis. Must be used with -list")

	var expectedNumArgs int
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

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

	if *matrixAverage != "" && *matrixComponents != "" {
		log.Fatalf("-matrixAverage and -matrixComponents cannot be used together")
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
	bedSimilarity(s)
}
