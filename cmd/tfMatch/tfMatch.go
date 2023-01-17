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
	OutFile            string
	MatrixFileType     string
	PropMatch          float64
	Pseudocounts       float64
	OutputasProportion bool
}

func tfMatch(s Settings) {
	records := fasta.Read(s.InFile)
	motifs := motif.ReadJaspar(s.MatrixFile, s.MatrixFileType)
	switch s.MatrixFileType {
	case "Frequency":
		motifs = motif.ReadJaspar(s.MatrixFile, "Frequency")
		motifs = motif.PfmSliceToPpmSlice(motifs, s.Pseudocounts)
		motifs = motif.PpmSliceToPwmSlice(motifs)
	case "Probability":
		motifs = motif.ReadJaspar(s.MatrixFile, "Probability")
		motifs = motif.PpmSliceToPwmSlice(motifs)
	case "Weight":
		motifs = motif.ReadJaspar(s.MatrixFile, "Weight")
	default:
		log.Fatalf("Error. Unexpected motif file format. Options are 'Frequency', 'Probability', and 'Weight'.")
	}
	if s.PropMatch < 0 || s.PropMatch > 1 {
		log.Fatalf("Error. PropMatch option should be a proportion, a value between 0 and 1.")
	}
	motif.RapidMatch(motifs, records, s.PropMatch, s.OutFile, s.OutputasProportion)
}

func usage() {
	fmt.Print(
		"tfMatch - Genome-wide scanning of TFBS occurrences." +
			"Input DNA sequences must be upper case.\n" +
			"Usage:\n" +
			"tfMatch input.fa matrices.pfm output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var propMatch *float64 = flag.Float64("propMatch", 0.8, "Specifies the minimum motif score (as a proportion of the consensus sequence score) required for a match to be retained in the output.")
	var matrixFileType *string = flag.String("matrixFileType", "Frequency", "Specify the type of position matrix file. Can be 'Frequency', 'Probability', or 'Weight'.")
	var pfmPseudocounts *float64 = flag.Float64("pfmPseudocounts", 0.1, "If a Position Frequency Matrix is provided, this pseudocount value will be applied when converting to a PWM.")
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
	outFile := flag.Arg(2)

	s := Settings{
		InFile:             inFile,
		MatrixFile:         matrixFile,
		OutFile:            outFile,
		MatrixFileType:     *matrixFileType,
		PropMatch:          *propMatch,
		Pseudocounts:       *pfmPseudocounts,
		OutputasProportion: *outputAsProportion,
	}

	tfMatch(s)
}
