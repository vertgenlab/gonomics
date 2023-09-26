package main

//TODO : add check for self-overlap bed files
//TODO : Remove file-name path in matrix output
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
	aNameSlice := strings.Split(s.bedA, "/")
	aName := aNameSlice[len(aNameSlice)-1]
	bNameSlice := strings.Split(s.bedB, "/")
	bName := bNameSlice[len(bNameSlice)-1]
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
	var a, b []bed.Bed
	var intervalsA, intervalsB []interval.Interval
	var percA, percB, avg float64
	var outMatrix *fileio.EasyWriter
	var err error
	var j int
	var l matrixLine
	var aNameSlice, bNameSlice []string
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
			aNameSlice = strings.Split(i, "/")
			aName = aNameSlice[len(aNameSlice)-1]
			allFiles = append(allFiles, aName)
		}
		header := strings.Join(allFiles, "\t")
		fileio.WriteToFileHandle(outMatrix, header)
	}

	for i := range files {
		aNameSlice = strings.Split(files[i], "/")
		aName = aNameSlice[len(aNameSlice)-1]
		l.name = aName
		l.vals = []float64{}
		a = bed.Read(files[i])
		intervalsA = interval.BedSliceToIntervals(a)
		for j = range files {
			if files[i] == files[j] {
				if s.matrixComponents != "" || s.matrixAverage != "" {
					l.vals = append(l.vals, 1.0)
				}
				continue
			}
			b = bed.Read(files[j])
			bNameSlice = strings.Split(files[j], "/")
			bName = bNameSlice[len(bNameSlice)-1]
			intervalsB = interval.BedSliceToIntervals(b)
			percA, percB, avg = interval.IntervalSimilarity(intervalsA, intervalsB)
			if j > i {
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

func bedSimilarityComp(s settings) {
	if s.list == "" {
		compareTwo(s)
	} else {
		multipleComparisons(s)
	}
}

func usage() {
	fmt.Print("haqerBedComps -- Takes in 2 bed files or a list of bed files and give similarity statistics based on number of overlaps" +
		"for the input bed files.\n" +
		"If the -list option is used, all combinations of provided files (excluding the comparison with itself) will be performed. " +
		"The statistics are: the percentage of elements that have an overlap an element in the second bed file, " +
		"the percentage of elements in the second bed file that overlap an element in the first bed file, and " +
		"the average of the two overlap percentages (a metric of overall similarity between the two bed files)." +
		"When using the -list option, either matrix option can be used to create a heatmap of similarities scores for all bed files in the list\n\n" +
		"Usage: \n" +
		"haqerBedComps inBedA inBedB outFile\n" +
		"or\n" +
		"haqerBedComps -list list.txt [options] outFile\n\n")
	flag.PrintDefaults()
}

func main() {
	var list *string = flag.String("list", "", "Provide a list of bed files and perform an IntervalSimilarity test on all possible combinations")
	var matrixAverage *string = flag.String("matrixAverage", "", "Provide a filename for a matrix that will have all files from the -list option on both axes. "+
		"The matrix will be populated with IntervalSimilarity metric. Must be used with -list")
	var matrixComponents *string = flag.String("matrixComponents", "", "Provide a filename for a matrix that will have all files from the -list option on both axes. "+
		"The matrix will be populated with the overlap percentage of elements from the file on the vertical axis with the file on the horizontal axis. Must be used with -list")

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
	bedSimilarityComp(s)
}
