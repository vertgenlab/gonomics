package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
	"log"
)

type InfoSettings struct {
	InFile       string
	OutFile      string
	MatrixType   string
	Pseudocounts float64
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
		records = motif.PfmSliceToPpmSlice(records, s.Pseudocounts)
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
