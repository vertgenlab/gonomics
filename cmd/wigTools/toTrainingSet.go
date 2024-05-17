package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"math"
	"math/rand"
	"os"
)

// ToTrainingSetSettings defines the usage settings for the wigTools toTrainingSet subcommand.
type ToTrainingSetSettings struct {
	InWigFile      string
	InFastaFile    string
	TrainFile      string
	ValidateFile   string
	TestFile       string
	WindowSize     int
	Stride         int
	ValidationProp float64
	TestingProp    float64
	SetSeed        int64
	Missing        float64
	LogTransform   bool
	IncludeRevComp bool
	NoHeader       bool
}

// toTrainingSetUsage defines the usage statement for the wigTools toTrainingSet subcommand.
func toTrainingSetUsage(toTrainingSetFlags *flag.FlagSet) {
	fmt.Printf("wigTools toTrainingSet - Converts a wig into training examples (sequence-amplitude pairs) for use in genomic sequence-to-function models.\n" +
		"This function can also partition these training examples into training, validation, and testing sets.\n" +
		"Usage:\n" +
		"wigTools toTrainingSet input.wig genome.fa train.txt validate.txt test.txt\n" +
		"options:\n")
	toTrainingSetFlags.PrintDefaults()
}

// parseToTrainingSetArgs is the main function of the wigTools toTrainingSet subcommand. It parses options and runs the toTrainingSet function.
func parseToTrainingSetArgs() {
	var expectedNumArgs int = 5
	var err error
	toTrainingSetFlags := flag.NewFlagSet("toTrainingSet", flag.ExitOnError)
	var missing *float64 = toTrainingSetFlags.Float64("missing", -10, "Sets the missing value. Regions with this number as the wig value will be skipped.")
	var windowSize *int = toTrainingSetFlags.Int("windowSize", 400, "Sets the windowSize for mapping sequence to scalars.")
	var stride *int = toTrainingSetFlags.Int("stride", 400, "Sets the stride between sequence sampling points. Recommended to be at least the windowSize to avoid correlation between training points.")
	var validationProp *float64 = toTrainingSetFlags.Float64("validationProp", 0.1, "Sets the proportion of the validation set. pValidate + pTesting + pTraining = 1.")
	var testingProp *float64 = toTrainingSetFlags.Float64("testingProp", 0.1, "Sets the proportion of the testing set. pValidate + pTesting + pTraining = 1.")
	var setSeed *int64 = toTrainingSetFlags.Int64("setSeed", -1, "Sets a seed for the random number generator.")
	var logTransform *bool = toTrainingSetFlags.Bool("logTransform", false, "Log transform the wig values in the output files.")
	var includeRevComp *bool = toTrainingSetFlags.Bool("includeRevComp", false, "Includes the rev comp for each sequence as an additional training/validation/testing example.")
	var noHeader *bool = toTrainingSetFlags.Bool("noHeader", false, "Do not include a header line in output data. Useful if subsetting or shuffling output.")
	err = toTrainingSetFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	toTrainingSetFlags.Usage = func() { toTrainingSetUsage(toTrainingSetFlags) }
	if len(toTrainingSetFlags.Args()) != expectedNumArgs {
		toTrainingSetFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d.\n", expectedNumArgs, len(toTrainingSetFlags.Args()))
	}

	inWigFile := toTrainingSetFlags.Arg(0)
	inFastaFile := toTrainingSetFlags.Arg(1)
	trainFile := toTrainingSetFlags.Arg(2)
	validateFile := toTrainingSetFlags.Arg(3)
	testFile := toTrainingSetFlags.Arg(4)

	s := ToTrainingSetSettings{
		InWigFile:      inWigFile,
		InFastaFile:    inFastaFile,
		TrainFile:      trainFile,
		ValidateFile:   validateFile,
		TestFile:       testFile,
		WindowSize:     *windowSize,
		Stride:         *stride,
		ValidationProp: *validationProp,
		TestingProp:    *testingProp,
		SetSeed:        *setSeed,
		Missing:        *missing,
		LogTransform:   *logTransform,
		IncludeRevComp: *includeRevComp,
		NoHeader:       *noHeader,
	}

	toTrainingSet(s)
}

// toTrainingSet creates training examples (sequence-amplitude pairs) from an input wig file and corresponding
// genome sequence. The output is split into a training, validation, and testing set.
func toTrainingSet(s ToTrainingSetSettings) {
	rand.Seed(s.SetSeed)
	var start, chromIndex, midpoint int
	var currName, lineToWrite string
	var currFa fasta.Fasta
	var currRand float64
	var err error

	trainOut := fileio.EasyCreate(s.TrainFile)
	testOut := fileio.EasyCreate(s.TestFile)
	validateOut := fileio.EasyCreate(s.ValidateFile)

	if !s.NoHeader {
		// write headers to all three output files
		_, err = fmt.Fprintf(trainOut, "name\tseq\tvalue\n")
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(testOut, "name\tseq\tvalue\n")
		exception.PanicOnErr(err)
		_, err = fmt.Fprintf(validateOut, "name\tseq\tvalue\n")
		exception.PanicOnErr(err)
	}

	if s.ValidationProp+s.TestingProp >= 1 {
		log.Fatalf("pVAlidation + pTesting should sum to less than one.")
	}

	wigs := wig.GoReadToChan(s.InWigFile)
	genome := fasta.Read(s.InFastaFile)

	for i := range wigs {
		chromIndex = getChromIndex(genome, i.Chrom)
		for start = 0; start < len(i.Values)-s.WindowSize; start += s.Stride {
			midpoint = (start + start + s.WindowSize) / 2
			if i.Values[midpoint] == s.Missing {
				continue //skip regions where there is no data in the wig.
			}
			currName = fmt.Sprintf("%s:%d-%d", i.Chrom, start, start+s.WindowSize)
			currFa = fasta.Extract(genome[chromIndex], start, start+s.WindowSize, currName)
			dna.AllToUpper(currFa.Seq)
			if s.LogTransform {
				lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", currFa.Name, dna.BasesToString(currFa.Seq), math.Log(i.Values[midpoint]))
			} else {
				lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", currFa.Name, dna.BasesToString(currFa.Seq), i.Values[midpoint]) //start+window / 2 is the midpoint of the window.
			}
			//now we shard the training example into either the testing, training, or validation set.
			currRand = rand.Float64()
			if currRand < s.TestingProp {
				_, err = fmt.Fprintf(testOut, lineToWrite)
				exception.PanicOnErr(err)
			} else if currRand < s.TestingProp+s.ValidationProp {
				_, err = fmt.Fprintf(validateOut, lineToWrite)
				exception.PanicOnErr(err)
			} else {
				_, err = fmt.Fprintf(trainOut, lineToWrite)
				exception.PanicOnErr(err)
			}

			if s.IncludeRevComp {
				dna.ReverseComplement(currFa.Seq)
				if s.LogTransform {
					lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", fmt.Sprintf("%s_rev", currFa.Name), dna.BasesToString(currFa.Seq), math.Log(i.Values[midpoint]))
				} else {
					lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", fmt.Sprintf("%s_rev", currFa.Name), dna.BasesToString(currFa.Seq), i.Values[midpoint])
				}
				//we don't regenerate currRand, for and rev for same sequence both land in the same shard.
				if currRand < s.TestingProp {
					_, err = fmt.Fprintf(testOut, lineToWrite)
					exception.PanicOnErr(err)
				} else if currRand < s.TestingProp+s.ValidationProp {
					_, err = fmt.Fprintf(validateOut, lineToWrite)
					exception.PanicOnErr(err)
				} else {
					_, err = fmt.Fprintf(trainOut, lineToWrite)
					exception.PanicOnErr(err)
				}
			}
		}
	}
	err = trainOut.Close()
	exception.PanicOnErr(err)
	err = testOut.Close()
	exception.PanicOnErr(err)
	err = validateOut.Close()
	exception.PanicOnErr(err)
}

// getChromIndex is a small helper function that finds the index in a []fasta.Fasta that
// has a queried Name, specified by the input string 'chrom'.
func getChromIndex(genome []fasta.Fasta, chrom string) int {
	for i := range genome {
		if genome[i].Name == chrom {
			return i
		}
	}
	log.Fatalf("Wig chromosome name: %s not found in target genome.", chrom)
	return -1
}
