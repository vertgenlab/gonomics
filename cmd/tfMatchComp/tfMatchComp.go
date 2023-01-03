package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/motif"
	"log"
)

type Settings struct {
	InFile             string
	MatrixFile         string
	ChromName          string
	OutFile            string
	PropMatch          float64
	MatrixFileType     string
	Pseudocounts       float64
	RefStart           int
	OutputAsProportion bool
}

func tfMatchComp(s Settings) {
	//read in fasta sequences
	records := fasta.Read(s.InFile)
	fasta.AllToUpper(records)

	// read and initialize position weight matrix
	var motifs []motif.PositionMatrix
	switch s.MatrixFileType {
	case "Frequency":
		motifs = motif.Read(s.MatrixFile, "Frequency")
		motifs = motif.PfmSliceToPpmSlice(motifs, s.Pseudocounts)
		motifs = motif.PpmSliceToPwmSlice(motifs)
	case "Probability":
		motifs = motif.Read(s.MatrixFile, "Probability")
		motifs = motif.PpmSliceToPwmSlice(motifs)
	case "Weight":
		motifs = motif.Read(s.MatrixFile, "Weight")
	default:
		log.Fatalf("Error. Unexpected motif file format. Options are 'Frequency', 'Probability', and 'Weight'.")
	}

	//pre-flight error checks
	if s.PropMatch < 0 || s.PropMatch > 1 {
		log.Fatalf("Error. PropMatch option should be a proportion, a value between 0 and 1.")
	}
	if len(records) != 2 {
		log.Fatalf("Error. tfMatchComp expects a pairwise multiFa alignment with two sequences. Found %v.\n", len(records))
	}
	if len(records[0].Seq) != len(records[1].Seq) {
		log.Fatalf("Error. tfMatchComp expects a well-formed pairwise multiFa alignment. Input sequences are not the same length. RefLen: %v. AltLen: %v.\n", len(records[0].Seq), len(records[1].Seq))
	}

	//run the program from the motif package
	motif.MatchComp(motifs, records, s.ChromName, s.PropMatch, s.OutFile, s.RefStart, s.OutputAsProportion)
}

func usage() {
	fmt.Print(
		"tfMatchComp - Compare the TFBS profiles between two input aligned genomic sequences.\n" +
			"Usage:\n" +
			"tfMatchComp input.fa matrices.pfm chromName output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var propMatch *float64 = flag.Float64("propMatch", 0.8, "Specifies the minimum motif score (as a proportion of the consensus sequence score) required for a match to be retained in the output.")
	var matrixFileType *string = flag.String("matrixFileType", "Frequency", "Specify the type of position matrix file. Can be 'Frequency', 'Probability', or 'Weight'.")
	var pfmPseudocounts *float64 = flag.Float64("pfmPseudocounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM.")
	var refStart *int = flag.Int("refStart", 0, "Set the reference position for the beginning of the input multiFa alignment.")
	var outputAsProportion *bool = flag.Bool("outputAsProportion", false, "Display the output motif scores as proportions of the consensus score. Motif difference score will thus be a change in consensus proportion.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	matrixFile := flag.Arg(1)
	chromName := flag.Arg(2)
	outFile := flag.Arg(3)

	s := Settings{
		InFile:             inFile,
		MatrixFile:         matrixFile,
		ChromName:          chromName,
		OutFile:            outFile,
		PropMatch:          *propMatch,
		MatrixFileType:     *matrixFileType,
		Pseudocounts:       *pfmPseudocounts,
		RefStart:           *refStart,
		OutputAsProportion: *outputAsProportion,
	}

	tfMatchComp(s)
}
