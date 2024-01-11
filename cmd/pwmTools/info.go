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

// InfoSettings defines the usage settings for the pwmTools info subcommand.
type InfoSettings struct {
	InFile       string
	OutFile      string
	MatrixType   string
	PseudoCounts float64
	GcContent    float64
	Threshold    float64
}

// infoUsage defines the usage statement for the pwmTools info subcommand.
func infoUsage(infoFlags *flag.FlagSet) {
	fmt.Printf("pwmTools info - a tool for reporting tabular statistics on PWM files.\n" +
		"Usage:\n" +
		"pwmTools info in.pwm out.txt\n" +
		"options:\n")
	infoFlags.PrintDefaults()
}

// parseInfoArgs is the main function of the pwmTools info subcommand. It parses options and initializes the pwmInfo function.
func parseInfoArgs() {
	var expectedNumArgs int = 2
	infoFlags := flag.NewFlagSet("info", flag.ExitOnError)
	var matrixType *string = infoFlags.String("matrixType", "Weight", "Specify the type of position matrix. May be one of: Frequency, Probability, or Weight.")
	var pfmPseudoCounts *float64 = infoFlags.Float64("pfmPseudoCounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM.")
	var gcContent *float64 = flag.Float64("gcContent", 0.5, "Set the expected GC content of the target sequence.")
	var threshold *float64 = infoFlags.Float64("threshold", 0.8, "Set the threshold value for motif matches. Motifs with scores above this value will be considered a match. Used for calculating cache size.")
	var err error
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
		GcContent:    *gcContent,
		Threshold:    *threshold,
	}

	pwmInfo(s)
}

// pwmInfo parses an input PositionMatrix file, produces summary statistics on each input matrix, and writes
// this information to an output tabular file.
// Output fields are TF (the motif identifier), MotifName (the human-readable motif name), Length (the length of
// the position matrix), Consensus Score (the position site scoring matrix value for the motif consensus sequence),
// and the CacheLength (the number of k-mers above a user-specified motif score).
func pwmInfo(s InfoSettings) {
	var err error
	var consensusSequence fasta.Fasta
	var consensusScore float64
	var couldScoreConsensusSequence bool
	var records []motif.PositionMatrix
	var currCache map[uint64]float64

	if s.Threshold > 1 || s.Threshold < 0 {
		log.Fatalf("Error: Threshold must be a value between 0 and 1. Found: %v.\n", s.Threshold)
	}

	switch s.MatrixType {
	case "Frequency":
		records = motif.ReadJaspar(s.InFile, "Frequency")
		records = motif.PfmSliceToPpmSlice(records, s.PseudoCounts)
		records = motif.PpmSliceToPwmSlice(records, s.GcContent)
	case "Probability":
		records = motif.ReadJaspar(s.InFile, "Probability")
		records = motif.PpmSliceToPwmSlice(records, s.GcContent)
	case "Weight":
		records = motif.ReadJaspar(s.InFile, "Weight")
	default:
		log.Fatalf("Error. Unexpected motif file format. Options are 'Frequency', 'Probability', and 'Weight'.")
	}

	out := fileio.EasyCreate(s.OutFile)

	_, err = fmt.Fprintf(out, "TF\tMotifName\tLength\tConsensusScore\tCacheLength\n")
	exception.PanicOnErr(err)

	for i := range records {
		consensusSequence = motif.ConsensusSequence(records[i], false)
		consensusScore, _, couldScoreConsensusSequence = motif.ScoreWindow(records[i], consensusSequence.Seq, 0)
		if !couldScoreConsensusSequence {
			log.Fatalf("Error: could not score consensus sequence for motif: %v.\n", records[i].Id)
		}
		currCache = motif.BuildKmerHash(records[i], s.Threshold)
		_, err = fmt.Fprintf(out, "%s\t%s\t%v\t%e\t%v\n", records[i].Name, records[i].Id, len(records[i].Mat[0]), consensusScore, len(currCache))
		exception.PanicOnErr(err)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}
