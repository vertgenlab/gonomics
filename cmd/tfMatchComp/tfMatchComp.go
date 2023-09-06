package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/motif"
	"log"
)

func tfMatchComp(s motif.MatchCompSettings, fastaFileName string) {
	//read in fasta sequences
	records := fasta.Read(fastaFileName)
	fasta.AllToUpper(records)
	s.Records = records

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
	motif.MatchComp(s)
}

func usage() {
	fmt.Print(
		"tfMatchComp - Compare the motif profiles between two input aligned genomic sequences." +
			"Output lines are as follows:\n" +
			"CHR\tCHROMSTART\tCHROMEND\tMOTIF_NAME\t0\tMOTIF_STRAND\tREF_SCORE\tALT_SCORE\tRESIDUAL" +
			"Usage:\n" +
			"tfMatchComp input.fa matrices.pfm/ppm/pwm chromName output.bed\n" +
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
	var residualWindowSize *int = flag.Int("residualWindowSize", 5, "Set the number of offset bases to consider when searching for motif in other species.")
	var enforceStrandMatch *bool = flag.Bool("enforceStrand", false, "If species A has the motif CCC, and the orthologous species B has sequence GGG, this is considered a match by default (as the motif is still there, just in revComp.\n"+
		"This option enforces strand matching.")
	var residualFilter *float64 = flag.Float64("residualFilter", 0, "The difference in motif scores between the two sequences must be at least this value to be retained in the output.")

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

	s := motif.MatchCompSettings{
		MotifFile:          matrixFile,
		MotifType:          *matrixFileType,
		PropMatch:          *propMatch,
		ChromName:          chromName,
		OutFile:            outFile,
		Pseudocounts:       *pfmPseudocounts,
		ResidualWindowSize: *residualWindowSize,
		RefStart:           *refStart,
		OutputAsProportion: *outputAsProportion,
		EnforceStrandMatch: *enforceStrandMatch,
		ResidualFilter:     *residualFilter,
	}

	tfMatchComp(s, inFile)
}
