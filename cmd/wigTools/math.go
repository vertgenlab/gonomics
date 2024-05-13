package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
)

// MathSettings defines the usage settings for the wigTools math subcommand.
type MathSettings struct {
	InFile                 string
	OutFile                string
	AbsoluteError          string
	AbsolutePercentError   string
	BedMask                string
	ChromSizes             string
	ElementWiseAdd         string
	ElementWiseMax         string
	ElementWiseSubtract    string
	MaxValue               float64
	MinValue               float64
	MovingAverageSmoothing int
	Missing                float64
	MissingBed             bool
	Pearson                string
	SamplingFrequency      float64
	ScalarMultiply         float64
	ScalarDivide           float64
	SetSeed                int64
}

// mathUsage defines the usage statement for the wigTools math subcommand.
func mathUsage(mathFlags *flag.FlagSet) {
	fmt.Print(
		"wigTools math - perform mathematical operations on wig format data.\n" +
			"Mathematical operations must be performed as single operations.\n" +
			"Usage:\n" +
			"wigTools math in.wig chrom.sizes out.file\n" +
			"options:\n")
	mathFlags.PrintDefaults()
}

// parseMathArgs is the main function of the wigTools math subcommand. It parses options and runs the wigMath function.
func parseMathArgs() {
	var expectedNumArgs = 3
	var err error
	mathFlags := flag.NewFlagSet("math", flag.ExitOnError)
	var absoluteError *string = mathFlags.String("absoluteError", "", "Specify a second wig file to determine the absolute error (element-wise) from the first. Returns a wig file.")
	var absolutePercentError *string = mathFlags.String("absolutePercentError", "", "Specify a second wig file to determine the absolute percent error (element-wise) from the first. Returns a wig file.")
	var bedMask *string = mathFlags.String("bedMask", "", "Specify a bed file, and mask all wig regions overlapping these bed regions to missing.")
	var elementWiseAdd *string = mathFlags.String("elementWiseAdd", "", "Specify a second wig file to add (element-wise) from the first. Returns a wig file.")
	var elementWiseMax *string = mathFlags.String("elementWiseMax", "", "Specify a second wig file to calculate the element-wise max value between the two wigs. Returns a wig file.")
	var elementWiseSubtract *string = mathFlags.String("elementWiseSubtract", "", "Specify a second wig file to subtract (element-wise) from the first. Returns a wig file.")
	var maxValue *float64 = mathFlags.Float64("maxValue", math.MaxFloat64, "Mask values in the output wig as 'missing' if they are above this value.")
	var minValue *float64 = mathFlags.Float64("minValue", math.MaxFloat64*-1, "Mask values in the output wig as 'missing' if they are below this value.")
	var missing *float64 = mathFlags.Float64("missing", 0, "Specify the value associated with missing data in the wig. These values are excluded from summary statistic calculations like Pearson.")
	var missingBed *bool = mathFlags.Bool("missingBed", false, "Create a bed file as the output file of all contiguous regions of the input wig with the missing value.")
	var movingAverageSmoothing *int = mathFlags.Int("movingAverageSmoothing", 1, "Set to a number greater than 1 to perform moving average smoothing on input wig data. Returns a wig file.")
	var pearson *string = mathFlags.String("pearson", "", "Specify a second wig file to determine the Pearson Correlation Coefficient between this wig and the first wig. Returns a text file.")
	var samplingFrequency *float64 = mathFlags.Float64("sampleFrequency", 0.001, "When calculating the Pearson correlation coefficient between two wigs, set the sampling frequency, or the proportion of positions in the wigs that will be evaluated in the calculation.")
	var scalarDivide *float64 = mathFlags.Float64("scalarDivide", 1, "Divide all entries in the input wig by the user-specified value (cannot divide by zero).")
	var scalarMultiply *float64 = mathFlags.Float64("scalarMultiply", 1, "Multiply all entries in the input wig by the user-specified value.")
	var setSeed *int64 = mathFlags.Int64("setSeed", 1, "Set the seed for the random number generator.")

	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	err = mathFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	mathFlags.Usage = func() { mathUsage(mathFlags) }
	if len(mathFlags.Args()) != expectedNumArgs {
		mathFlags.Usage()
		log.Fatalf("Error: expected %d arguments, but got %d.\n", expectedNumArgs, len(mathFlags.Args()))
	}
	inFile := mathFlags.Arg(0)
	chromSizes := mathFlags.Arg(1)
	outFile := mathFlags.Arg(2)

	s := MathSettings{
		InFile:                 inFile,
		OutFile:                outFile,
		AbsoluteError:          *absoluteError,
		AbsolutePercentError:   *absolutePercentError,
		BedMask:                *bedMask,
		ChromSizes:             chromSizes,
		ElementWiseAdd:         *elementWiseAdd,
		ElementWiseMax:         *elementWiseMax,
		ElementWiseSubtract:    *elementWiseSubtract,
		MinValue:               *minValue,
		MaxValue:               *maxValue,
		MovingAverageSmoothing: *movingAverageSmoothing,
		Missing:                *missing,
		MissingBed:             *missingBed,
		Pearson:                *pearson,
		SamplingFrequency:      *samplingFrequency,
		ScalarMultiply:         *scalarMultiply,
		ScalarDivide:           *scalarDivide,
		SetSeed:                *setSeed,
	}

	wigMath(s)
}

// wigMth performs mathematical operations on input wig format files.
func wigMath(s MathSettings) {
	rand.Seed(s.SetSeed)
	var err error
	var second map[string]wig.Wig
	var foundInMap bool
	var currPos int
	var currKey string
	records := wig.Read(s.InFile, s.ChromSizes, s.Missing)
	multipleOptionCheck(s)
	if s.ScalarMultiply != 1 {
		for currKey = range records {
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] *= s.ScalarMultiply
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ScalarDivide != 1 {
		if s.ScalarDivide == 0 {
			log.Fatalf("You cannot divide wig values by zero, please input another value. I'm not mad but I'm a bit disappointed.")
		}
		for currKey = range records {
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] /= s.ScalarDivide
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ElementWiseAdd != "" {
		second = wig.Read(s.ElementWiseAdd, s.ChromSizes, s.Missing)
		for currKey = range records {
			if _, foundInMap = second[currKey]; !foundInMap {
				log.Fatalf("Error: Chrom in first wig file: %v, not found in the second wig file.\n", currKey)
			}
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && second[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] += second[currKey].Values[currPos]
				} else {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ElementWiseMax != "" {
		second = wig.Read(s.ElementWiseMax, s.ChromSizes, s.Missing)
		for currKey = range records {
			if _, foundInMap = second[currKey]; !foundInMap {
				log.Fatalf("Error: Chrom in first wig file: %v, not found in the second wig file.\n", currKey)
			}
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && second[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] = numbers.Max(records[currKey].Values[currPos], second[currKey].Values[currPos])
				} else {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.ElementWiseSubtract != "" {
		second = wig.Read(s.ElementWiseSubtract, s.ChromSizes, s.Missing)
		for currKey = range records {
			if _, foundInMap = second[currKey]; !foundInMap {
				log.Fatalf("Error: Chrom in first wig file: %v, not found in the second wig file.\n", currKey)
			}
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && second[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] -= second[currKey].Values[currPos]
				} else {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.MovingAverageSmoothing > 1 {
		records = wig.SmoothMap(records, s.MovingAverageSmoothing, s.Missing)
		wig.Write(s.OutFile, records)
	} else if s.AbsoluteError != "" {
		second = wig.Read(s.AbsoluteError, s.ChromSizes, s.Missing)
		for currKey = range records {
			if _, foundInMap = second[currKey]; !foundInMap {
				log.Fatalf("Error: Chrom in first wig file: %v, not found in the second wig file.\n", currKey)
			}
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && second[currKey].Values[currPos] != s.Missing {
					records[currKey].Values[currPos] = math.Abs(records[currKey].Values[currPos] - second[currKey].Values[currPos])
				} else {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.AbsolutePercentError != "" {
		second = wig.Read(s.AbsolutePercentError, s.ChromSizes, s.Missing)
		for currKey = range records {
			for currPos = range records[currKey].Values {
				if _, foundInMap = second[currKey]; !foundInMap {
					log.Fatalf("Error: Chrom in first wig file: %v, not found in the second wig file.\n", currKey)
				}
				if records[currKey].Values[currPos] == s.Missing || second[currKey].Values[currPos] == s.Missing {
					records[currKey].Values[currPos] = s.Missing
				} else if records[currKey].Values[currPos] != 0 {
					records[currKey].Values[currPos] = math.Abs((records[currKey].Values[currPos]-second[currKey].Values[currPos])/records[currKey].Values[currPos]) * 100
				} else {
					records[currKey].Values[currPos] = s.Missing //these positions are undefined.
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.Pearson != "" {
		second = wig.Read(s.Pearson, s.ChromSizes, s.Missing)
		answer := wig.Pearson(records, second, s.Missing, s.SamplingFrequency)
		out := fileio.EasyCreate(s.OutFile)
		_, err = fmt.Fprintf(out, "PCC:\t%f\n", answer)
		err = out.Close()
		exception.PanicOnErr(err)
	} else if s.MinValue > -1*math.MaxFloat64 {
		for currKey = range records {
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && records[currKey].Values[currPos] < s.MinValue {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.MaxValue < math.MaxFloat64 {
		for currKey = range records {
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] != s.Missing && records[currKey].Values[currPos] > s.MaxValue {
					records[currKey].Values[currPos] = s.Missing
				}
			}
		}
		wig.Write(s.OutFile, records)
	} else if s.MissingBed {
		var inMissingRegion = false
		out := fileio.EasyCreate(s.OutFile)
		var currentBed = bed.Bed{Chrom: "dummyPlaceHolder", ChromStart: -1, ChromEnd: -1, FieldsInitialized: 3}

		//to make this program deterministic, we read the keys to a slice, sort the slice, and then iterate through the map in sorted key order
		keys := make([]string, 0)
		for currKey = range records {
			keys = append(keys, currKey)
		}
		sort.Strings(keys)

		for _, currKey = range keys {
			for currPos = range records[currKey].Values {
				if records[currKey].Values[currPos] == s.Missing {
					if records[currKey].Chrom != currentBed.Chrom && currentBed.Chrom != "dummyPlaceHolder" { //case where previous chrom ended in missing,
						// so we write the bed entry, but new chrom also starts with missing.
						bed.WriteBed(out, currentBed)
						currentBed = bed.Bed{Chrom: records[currKey].Chrom, ChromStart: currPos, ChromEnd: currPos + 1, FieldsInitialized: 3}
					} else if inMissingRegion {
						currentBed.ChromEnd = currPos + 1 //include currentPos in currentBed.
					} else {
						currentBed = bed.Bed{Chrom: records[currKey].Chrom, ChromStart: currPos, ChromEnd: currPos + 1, FieldsInitialized: 3}
						inMissingRegion = true
					}
				} else {
					if inMissingRegion {
						inMissingRegion = false
						bed.WriteBed(out, currentBed)
					}
				}
			}
		}
		//if we found a valid bed, we write the last one out
		if currentBed.ChromStart >= 0 {
			bed.WriteBed(out, currentBed)
		}
		err = out.Close()
		exception.PanicOnErr(err)
	} else if s.BedMask != "" {
		var currBed int
		bedRegions := bed.Read(s.BedMask)
		for currBed = range bedRegions {
			for currPos = bedRegions[currBed].ChromStart; currPos < bedRegions[currBed].ChromEnd; currPos++ {
				if currPos >= len(records[bedRegions[currBed].Chrom].Values) {
					log.Fatalf("Error: current position (%v) exceeds length of chromosome %s.\n", currPos, bedRegions[currBed].Chrom)
				}
				records[bedRegions[currBed].Chrom].Values[currPos] = s.Missing
			}
		}
		wig.Write(s.OutFile, records)
	}
}

func multipleOptionCheck(s MathSettings) {
	var optionCount = 0
	if s.BedMask != "" {
		optionCount++
	}
	if s.MinValue > -1*math.MaxFloat64 {
		optionCount++
	}
	if s.MaxValue < math.MaxFloat64 {
		optionCount++
	}
	if s.ScalarMultiply != 1 {
		optionCount++
	}
	if s.ElementWiseAdd != "" {
		optionCount++
	}
	if s.ElementWiseMax != "" {
		optionCount++
	}
	if s.ElementWiseSubtract != "" {
		optionCount++
	}
	if s.MissingBed {
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
		log.Fatalf("wigTools math can perform only one mathematical operation at a time. Rerun with a single option.")
	}
}
