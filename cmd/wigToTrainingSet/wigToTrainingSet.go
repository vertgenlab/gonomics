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
)

func wigToTrainingSet(s Settings) {
	rand.Seed(s.SetSeed)
	var start, chromIndex int
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

	if s.ValidationProp+s.TestingProp > 1 {
		log.Fatalf("pVAlidation + pTesting should sum to less than one.")
	}

	wigs := wig.GoReadToChan(s.InWigFile)
	genome := fasta.Read(s.InFastaFile)

	for i := range wigs {
		chromIndex = getChromIndex(genome, i.Chrom)
		for start = 0; start < len(i.Values)-s.WindowSize; start += s.Stride {
			if i.Values[(start+start+s.WindowSize)/2] == s.Missing {
				continue //skip regions where there is no data in the wig.
			}
			currName = fmt.Sprintf("%s:%v-%v", i.Chrom, start, start+s.WindowSize)
			currFa = fasta.Extract(genome[chromIndex], start, start+s.WindowSize, currName)
			dna.AllToUpper(currFa.Seq)
			if s.LogTransform {
				lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", currFa.Name, dna.BasesToString(currFa.Seq), math.Log(i.Values[(start+start+s.WindowSize)/2]))
			} else {
				lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", currFa.Name, dna.BasesToString(currFa.Seq), i.Values[(start+start+s.WindowSize)/2]) //start+window / 2 is the midpoint of the window.
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
					lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", fmt.Sprintf("%s_rev", currFa.Name), dna.BasesToString(currFa.Seq), math.Log(i.Values[(start+start+s.WindowSize)/2]))
				} else {
					lineToWrite = fmt.Sprintf("%s\t%s\t%g\n", fmt.Sprintf("%s_rev", currFa.Name), dna.BasesToString(currFa.Seq), i.Values[(start+start+s.WindowSize)/2])
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

func getChromIndex(genome []fasta.Fasta, chrom string) int {
	for i := range genome {
		if genome[i].Name == chrom {
			return i
		}
	}
	log.Fatalf("Wig chromosome name: %v not found in target genome.", chrom)
	return -1
}

func usage() {
	fmt.Print(
		"wigToTrainingSet - Converts a wig to a training, validation, and testing set for the GenomeSequenceConvNet.\n" +
			"Usage:\n" +
			"wigToTrainingSet input.wig genome.fa train.txt validate.txt test.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
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
	NoHeader bool
}

func main() {
	var expectedNumArgs int = 5
	var missing *float64 = flag.Float64("missing", -10, "Sets the missing value. Regions with this number as the wig value will be skipped.")
	var windowSize *int = flag.Int("windowSize", 400, "Sets the windowSize for mapping sequence to scalars.")
	var stride *int = flag.Int("stride", 400, "Sets the stride between sequence sampling points. Recommended to be at least the windowSize to avoid correlation between training points.")
	var validationProp *float64 = flag.Float64("validationProp", 0.1, "Sets the proportion of the validation set. pValidate + pTesting + pTraining = 1.")
	var testingProp *float64 = flag.Float64("testingProp", 0.1, "Sets the proportion of the testing set. pValidate + pTesting + pTraining = 1.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Sets a seed for the random number generator.")
	var logTransform *bool = flag.Bool("logTransform", false, "Log transform the wig values in the output files.")
	var includeRevComp *bool = flag.Bool("includeRevComp", false, "Includes the rev comp for each sequence as an additional training/validation/testing example.")
	var noHeader *bool = flag.Bool("noHeader", false, "Do not include a header line in output data. Useful if subsetting or shuffling output.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inWigFile := flag.Arg(0)
	inFastaFile := flag.Arg(1)
	trainFile := flag.Arg(2)
	validateFile := flag.Arg(3)
	testFile := flag.Arg(4)

	s := Settings{
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
		NoHeader: *noHeader,
	}

	wigToTrainingSet(s)
}
