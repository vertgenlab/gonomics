package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"math/rand"
)

func wigMath(s Settings) {
	rand.Seed(s.SetSeed)
	var err error
	var second []wig.Wig
	var i, j, chromIndex int
	records := wig.Read(s.InFile)
	multipleOptionCheck(s)
	if s.ScalarMultiply != 1 {
		for i = range records {
			for j = range records[i].Values {
				records[i].Values[j] *= s.ScalarMultiply
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ScalarDivide != 1 {
		if s.ScalarDivide == 0 {
			log.Fatalf("You cannot divide wig values by zero, please input another value. I'm not mad but I'm a bit disappointed.")
		}
		for i = range records {
			for j = range records[i].Values {
				records[i].Values[j] /= s.ScalarDivide
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ElementWiseAdd != "" {
		second = wig.Read(s.ElementWiseAdd)
		for i = range records {
			chromIndex = getChromIndex(second, records[i].Chrom)
			for j = range records[i].Values {
				records[i].Values[j] += second[chromIndex].Values[j]
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ElementWiseSubtract != "" {
		second = wig.Read(s.ElementWiseSubtract)
		for i = range records {
			chromIndex = getChromIndex(second, records[i].Chrom)
			for j = range records[i].Values {
				records[i].Values[j] -= second[chromIndex].Values[j]
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.MovingAverageSmoothing > 1 {
		records = wig.SmoothSlice(records, s.MovingAverageSmoothing)
		wig.Write(s.OutFile, records)
	} else if s.AbsoluteError != "" {
		second = wig.Read(s.AbsoluteError)
		for i = range records {
			chromIndex = getChromIndex(second, records[i].Chrom)
			for j = range records[i].Values {
				records[i].Values[j] = math.Abs(records[i].Values[j] - second[chromIndex].Values[j])
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.AbsolutePercentError != "" {
		second = wig.Read(s.AbsolutePercentError)
		for i = range records {
			chromIndex = getChromIndex(second, records[i].Chrom)
			for j = range records[i].Values {
				if records[i].Values[j] != 0 {
					records[i].Values[j] = math.Abs((records[i].Values[j]-second[chromIndex].Values[j])/records[i].Values[j]) * 100
				} else {
					records[i].Values[j] = 0 //placeholder, these positions are undefined.
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.Pearson != "" {
		second = wig.Read(s.Pearson)
		answer := wig.Pearson(records, second, s.Missing, s.SamplingFrequency)
		out := fileio.EasyCreate(s.OutFile)
		_, err = fmt.Fprintf(out, "PCC:\t%f\n", answer)
		err = out.Close()
		exception.PanicOnErr(err)
	}
}

func multipleOptionCheck(s Settings) {
	var optionCount = 0
	if s.ScalarMultiply != 1 {
		optionCount++
	}
	if s.ElementWiseAdd != "" {
		optionCount++
	}
	if s.ElementWiseSubtract != "" {
		optionCount++
	}
	if s.MovingAverageSmoothing > 1 {
		optionCount++
	}
	if s.AbsoluteError != "" {
		optionCount++
	}
	if s.AbsolutePercentError != "" {
		optionCount++
	}
	if s.Pearson != "" {
		optionCount++
	}
	if optionCount > 1 {
		log.Fatalf("wigMath can perform only one mathematical operation at a time. Rerun with a single option.")
	}
}

func getChromIndex(w []wig.Wig, chrom string) int {
	for i := range w {
		if w[i].Chrom == chrom {
			return i
		}
	}
	log.Fatalf("Error. Chromosome %v not found in wig.", chrom)
	return -1
}

func usage() {
	fmt.Print(
		"wigMath - Perform mathematical operations on wig format data.\n" +
			"Mathematical operations must be performed as single operations.\n" +
			"Usage:\n" +
			"wigMath in.wig out.file" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile                 string
	OutFile                string
	ScalarMultiply float64
	ScalarDivide float64
	ElementWiseAdd string
	ElementWiseSubtract    string
	MovingAverageSmoothing int
	AbsoluteError          string
	AbsolutePercentError   string
	Missing                float64
	Pearson                string
	SamplingFrequency      float64
	SetSeed                int64
}

func main() {
	var expectedNumArgs = 2
	var scalarMultiply *float64 = flag.Float64("scalarMultiply", 1, "Multiply all entries in the input wig by the user-specified value.")
	var scalarDivide *float64 = flag.Float64("scalarDivide", 1, "Divide all entries in the input wig by the user-specified value (cannot divide by zero).")
	var elementWiseAdd *string = flag.String("elementWiseAdd", "", "Specify a second wig file to add (element-wise) from the first. Returns a wig file.")
	var elementWiseSubtract *string = flag.String("elementWiseSubtract", "", "Specify a second wig file to subtract (element-wise) from the first. Returns a wig file.")
	var movingAverageSmoothing *int = flag.Int("movingAverageSmoothing", 1, "Set to a number greater than 1 to perform moving average smoothing on input wig data. Returns a wig file.")
	var absoluteError *string = flag.String("absoluteError", "", "Specify a second wig file to determine the absolute error (element-wise) from the first. Returns a wig file.")
	var absolutePercentError *string = flag.String("absolutePercentError", "", "Specify a second wig file to determine the absolute percent error (element-wise) from the first. Returns a wig file.")
	var missing *float64 = flag.Float64("missing", 0, "Specify the value associated with missing data in the wig. These values are excluded from summary statistic calculations like Pearson.")
	var pearson *string = flag.String("pearson", "", "Specify a second wig file to determine the Pearson Correlation Coefficient between this wig and the first wig. Returns a text file.")
	var samplingFrequency *float64 = flag.Float64("sampleFrequency", 0.001, "When calculating the Pearson correlation coefficient between two wigs, set the sampling frequency, or the proportion of positions in the wigs that will be evaluated in the calculation.")
	var setSeed *int64 = flag.Int64("setSeed", 1, "Set the seed for the random number generator.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)
	s := Settings{
		InFile:                 inFile,
		OutFile:                outFile,
		ScalarMultiply: *scalarMultiply,
		ScalarDivide: *scalarDivide,
		ElementWiseAdd: *elementWiseAdd,
		ElementWiseSubtract:    *elementWiseSubtract,
		MovingAverageSmoothing: *movingAverageSmoothing,
		AbsoluteError:          *absoluteError,
		AbsolutePercentError:   *absolutePercentError,
		Missing:                *missing,
		Pearson:                *pearson,
		SamplingFrequency:      *samplingFrequency,
		SetSeed:                *setSeed,
	}
	wigMath(s)
}
