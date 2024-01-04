package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
	"log"
	"math/rand"
	"os"
)

// ShuffleSettings defines the usage options for the pwmTools shuffle subcommand.
type ShuffleSettings struct {
	InFile     string
	OutFile    string
	NumShuffle int
	SetSeed    int64
}

// shuffleUsage defines the usage statement for the pwmTools shuffle subcommand.
func shuffleUsage(shuffleFlags *flag.FlagSet) {
	fmt.Printf("pwmTools shuffle - a tool for producing shuffled control motifs from input motifs.\n" +
		"Input motifs are shuffled column-wise, preserving the base composition but randomizing motif order.\n" +
		"Useful for producing control motifs for excess conservation analysis.\n" +
		"pwmTools shuffle in.pwm out.pwm\n" +
		"options:\n")
	shuffleFlags.PrintDefaults()
}

// this is equivalent to the main function for the shuffle subcommand.
func parseShuffleArgs() {
	var expectedNumArgs int = 2
	var err error
	shuffleFlags := flag.NewFlagSet("shuffle", flag.ExitOnError)
	var numShuffle *int = shuffleFlags.Int("numShuffle", 10, "Specify the number of desired shuffled motifs per input motif.")
	var setSeed *int64 = shuffleFlags.Int64("setSeed", 1, "Specify the seed for the random number generator.")
	err = shuffleFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	shuffleFlags.Usage = func() { shuffleUsage(shuffleFlags) }

	if len(shuffleFlags.Args()) != expectedNumArgs {
		shuffleFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(shuffleFlags.Args()))
	}

	inFile := shuffleFlags.Arg(0)
	outFile := shuffleFlags.Arg(1)

	s := ShuffleSettings{
		InFile:     inFile,
		OutFile:    outFile,
		NumShuffle: *numShuffle,
		SetSeed:    *setSeed,
	}

	pwmShuffle(s)
}

// pwmShuffle is the primary function for this subcommand. It first parses an input file of PositionMatrices.
// For each input motif, a number of shuffled motifs are generated (default 10) and then written to an output file.
func pwmShuffle(s ShuffleSettings) {
	var originalMotifName string
	var currIter int
	var err error
	rand.Seed(s.SetSeed)
	records := motif.ReadJaspar(s.InFile, "Frequency") // we don't use the type information, so we'll just read it as Frequency, though it doesn't matter
	out := fileio.EasyCreate(s.OutFile)
	for currMatrix := range records {
		originalMotifName = records[currMatrix].Name
		for currIter = 0; currIter < s.NumShuffle; currIter++ {
			shufflePwmColumns(records[currMatrix])
			records[currMatrix].Name = fmt.Sprintf("%s_%v", originalMotifName, currIter)
			motif.WritePositionMatrixJaspar(out, records[currMatrix])
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

// shufflePwmColumns randomizes the order of columns in an input PositionMatrix in place using the Fisher-Yates algorithm.
func shufflePwmColumns(pwm motif.PositionMatrix) {
	var currRand, currRow int
	for currCol := range pwm.Mat[0] {
		currRand = rand.Intn(currCol + 1)
		for currRow = range pwm.Mat {
			pwm.Mat[currRow][currCol], pwm.Mat[currRow][currRand] = pwm.Mat[currRow][currRand], pwm.Mat[currRow][currCol]
		}
	}
}
