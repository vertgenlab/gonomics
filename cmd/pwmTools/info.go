package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
	"log"
	"os"
)

type InfoSettings struct {
	InFile       string
	OutFile      string
	MatrixType   string
	PseudoCounts float64
}

func infoUsage(infoFlags *flag.FlagSet) {
	fmt.Printf("pwmTools info - a tool for reporting tabular statistics on PWM files.\n" +
		"Usage:\n" +
		"pwmTools info in.pwm out.txt\n" +
		"options:\n")
	infoFlags.PrintDefaults()
}

func parseInfoArgs() {
	var err error
	var expectedNumArgs int = 2
	infoFlags := flag.NewFlagSet("info", flag.ExitOnError)

	var matrixType *string = infoFlags.String("matrixType", "Weight", "Specify the type of position matrix. May be one of: Frequency, Probability, or Weight.")
	var pfmPseudoCounts *float64 = infoFlags.Float64("pfmPseudoCounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM.")
	err = infoFlags.Parse(os.Args[2:])
	exception.PanicOnErr(err)
	infoFlags.Usage = func() { infoUsage(infoFlags) }
	if len(infoFlags.Args()) != expectedNumArgs {
		infoFlags.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(infoFlags.Args()))
	}

	inFile := infoFlags.Arg(0)
	outFile := infoFlags.Arg(1)

	s := InfoSettings{
		InFile:       inFile,
		OutFile:      outFile,
		MatrixType:   *matrixType,
		PseudoCounts: *pfmPseudoCounts,
	}

	pwmInfo(s)
}

func pwmInfo(s InfoSettings) {
	var err error
	var consensusSequence fasta.Fasta
	var consensusScore float64
	var couldScoreConsensusSequence bool
	var records []motif.PositionMatrix
	switch s.MatrixType {
	case "Frequency":
		records = motif.ReadJaspar(s.InFile, "Frequency")
		records = motif.PfmSliceToPpmSlice(records, s.PseudoCounts)
		records = motif.PpmSliceToPwmSlice(records)
	case "Probability":
		records = motif.ReadJaspar(s.InFile, "Probability")
		records = motif.PpmSliceToPwmSlice(records)
	case "Weight":
		records = motif.ReadJaspar(s.InFile, "Weight")
	default:
		log.Fatalf("Error. Unexpected motif file format. Options are 'Frequency', 'Probability', and 'Weight'.")
	}

	out := fileio.EasyCreate(s.OutFile)

	_, err = fmt.Fprintf(out, "TF\tMotifName\tLength\tConsensusScore\n")
	exception.PanicOnErr(err)

	for i := range records {
		consensusSequence = motif.ConsensusSequence(records[i], false)
		consensusScore, _, couldScoreConsensusSequence = motif.ScoreWindow(records[i], consensusSequence.Seq, 0)
		if !couldScoreConsensusSequence {
			log.Fatalf("Error: could not score consensus sequence for motif: %v.\n", records[i].Id)
		}
		_, err = fmt.Fprintf(out, "%s\t%s\t%v\t%e\n", records[i].Name, records[i].Id, len(records[i].Mat[0]), consensusScore)
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
