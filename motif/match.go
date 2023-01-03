package motif

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

// MatchComp computes all TFBS matches in a pairwise multiFa alignment for all motifs defined in a []PositionMatrix
// above an input propMatch (the motif score measured as a proportion of the consensus motif score).
// MatchComp calculates the motif score for both the ref and alt sequence, as well as the difference in motif score
// between the two sequences, and reports these values for all identified motifs in the additional fields of output bed
// entries.
func MatchComp(motifs []PositionMatrix, records []fasta.Fasta, chromName string, propMatch float64, outFile string, refStart int, outputAsProportion bool) {
	var currConsensus fasta.Fasta
	var currConsensusScore, currRefScore, currAltScore, currMotifScoreDifference float64
	var couldScoreRefWindow, couldScoreAltWindow, couldScoreConsensus bool
	var currCutoff float64
	var currAnnotation []string
	var currOutput bed.Bed
	var currReverseMotif PositionMatrix
	var err error

	out := fileio.EasyCreate(outFile)

	for _, motif := range motifs {
		var refPos int = refStart
		currConsensus = ConsensusSequence(motif, false)                                    //tieBreak can be false because either way we'll have the same score.
		currConsensusScore, couldScoreConsensus = ScoreWindow(motif, currConsensus.Seq, 0) //consensus sequence is necessarily positive strand by convention.
		if !couldScoreConsensus {
			log.Fatalf("Error scoring consensus sequence.")
		}
		currCutoff = propMatch * currConsensusScore
		currReverseMotif = ReverseComplement(motif)

		for alnPos := range records[0].Seq {
			if dna.DefineBase(records[0].Seq[alnPos]) {
				currRefScore, couldScoreRefWindow = ScoreWindow(motif, records[0].Seq, alnPos)
				currAltScore, couldScoreAltWindow = ScoreWindow(motif, records[1].Seq, alnPos)
				if couldScoreRefWindow && couldScoreAltWindow && (currRefScore > currCutoff || currAltScore > currCutoff) {
					if outputAsProportion {
						currMotifScoreDifference = (currRefScore / currConsensusScore) - (currAltScore / currConsensusScore)
						currAnnotation = []string{fmt.Sprintf("%f", currRefScore/currConsensusScore), fmt.Sprintf("%f", currAltScore/currConsensusScore), fmt.Sprintf("%f", currMotifScoreDifference)}
					} else {
						currMotifScoreDifference = currRefScore - currAltScore
						currAnnotation = []string{fmt.Sprintf("%f", currRefScore), fmt.Sprintf("%f", currAltScore), fmt.Sprintf("%f", currMotifScoreDifference)}
					}
					currOutput = bed.Bed{Chrom: chromName,
						ChromStart:        refPos,
						ChromEnd:          refPos + len(motif.Mat[0]),
						Name:              motif.Name,
						Score:             0,
						Strand:            bed.Positive,
						FieldsInitialized: 9,
						Annotation:        currAnnotation}
					bed.WriteBed(out, currOutput)
				}

				// now we score the reverse strand
				currRefScore, couldScoreRefWindow = ScoreWindow(currReverseMotif, records[0].Seq, alnPos)
				currAltScore, couldScoreAltWindow = ScoreWindow(currReverseMotif, records[1].Seq, alnPos)
				if couldScoreRefWindow && couldScoreAltWindow && (currRefScore > currCutoff || currAltScore > currCutoff) {
					if outputAsProportion {
						currMotifScoreDifference = (currRefScore / currConsensusScore) - (currAltScore / currConsensusScore)
						currAnnotation = []string{fmt.Sprintf("%f", currRefScore/currConsensusScore), fmt.Sprintf("%f", currAltScore/currConsensusScore), fmt.Sprintf("%f", currMotifScoreDifference)}
					} else {
						currMotifScoreDifference = currRefScore - currAltScore
						currAnnotation = []string{fmt.Sprintf("%f", currRefScore), fmt.Sprintf("%f", currAltScore), fmt.Sprintf("%f", currMotifScoreDifference)}
					}
					currOutput = bed.Bed{Chrom: chromName,
						ChromStart:        refPos,
						ChromEnd:          refPos + len(currReverseMotif.Mat[0]),
						Name:              currReverseMotif.Name,
						Score:             0,
						Strand:            bed.Negative,
						FieldsInitialized: 9,
						Annotation:        currAnnotation}
					bed.WriteBed(out, currOutput)
				}
				refPos++
			}
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

// bool returns false if the sequence could not be scored (contains Ns, ran off end of sequence)
func ScoreWindow(pm PositionMatrix, seq []dna.Base, alnStart int) (float64, bool) {
	var currAlnPos = alnStart
	var answer float64
	var motifPos int = 0
	for motifPos < len(pm.Mat[0]) {
		if currAlnPos >= len(seq) {
			return -1, false //ran off the end of the sequence.
		}
		switch seq[currAlnPos] {
		case dna.Gap:
			//we advance in alnPos at the end of the loop to check the next position but do not advance the motifPosition.
		case dna.A:
			answer += pm.Mat[0][motifPos]
			motifPos++
		case dna.C:
			answer += pm.Mat[1][motifPos]
			motifPos++
		case dna.G:
			answer += pm.Mat[2][motifPos]
			motifPos++
		case dna.T:
			answer += pm.Mat[3][motifPos]
			motifPos++
		case dna.N:
			return -1, false //this will just ignore this position of the sequence and move on.
		default:
			log.Fatalf("Unrecognized base. Cannot score window.")
			return -1, false
		}
		currAlnPos++
	}
	return answer, true
}
