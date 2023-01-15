package motif

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaTwoBit"
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

// RankTensorElement encodes the value and base for a RankTensor at a particular row and column.
// In three dimensions, it is useful to consider the Value and Base as two layers of the RankTensor.
type RankTensorElement struct {
	Value float64
	Base dna.Base
}

func initializeRankTensor(p PositionMatrix) [][]RankTensorElement {
	var row, column int
	var currMaxRow, currRank int
	var currMaxValue float64
	var answer = make([][]RankTensorElement, 4)
	for row = 0; row < 4; row ++ {
		answer[row] = make([]RankTensorElement, len(p.Mat[row]))
		for column = 0; column < len(p.Mat[row]); column++ {
			answer[row][column] = RankTensorElement{Value: p.Mat[row][column], Base: dna.Base(row)}//cast row to dna.Base (0 -> A, 1 -> C, 2 -> G, 3 -> T)
		}
	}
	//sort matrix columns by rank
	for column = 0; column < len(p.Mat[0]); column++ {
		for currRank = 0; currRank < 3; currRank++ {
			currMaxRow = currRank
			currMaxValue = answer[currRank][column].Value
			for row = currRank+1; row < 4; row++ {
				if answer[row][column].Value > currMaxValue {
					currMaxRow = row
					currMaxValue = answer[row][column].Value
				}
			}
			//swap max value with value in current rank row
			answer[currMaxRow][column], answer[currRank][column] = answer[currRank][column], answer[currMaxRow][column]
		}
	}
	return answer
}

// rankTensorToString formats a RankTensor as a string for debugging and visualization.
func rankTensorToString(m [][]RankTensorElement) string {
	var row, column int
	var answer,currBaseString string
	for row = 0; row < len(m); row++ {
		answer = answer + "[\t"
		for column = 0; column < len(m[row]); column++ {
			currBaseString = dna.BaseToString(m[row][column].Base)
			answer = answer + fmt.Sprintf("%f:%s\t", m[row][column].Value, currBaseString)
		}
		answer = answer + "]\n"
	}
	answer = answer + "\n"

	return answer
}

// buildKmerHas produces a hash mapping 2bit encoded kmer sequences to their corresponding motif score for a PositionMatrix.
// Only kmers with a motif score above an input thresholdProportion are stored in the
func buildKmerHash(p PositionMatrix, thresholdProportion float64) map[uint64]float64 {
	var answer = make(map[uint64]float64)
	var currSeq = ConsensusSequence(p, false)
	consensusValue, couldScoreConsensus := ScoreWindow(p, currSeq.Seq, 0)
	if !couldScoreConsensus {
		log.Fatalf("Error in buildKmerHash. Could not score consensus sequence.")
	}
	var threshold float64 = thresholdProportion * consensusValue
	var rankMatrix [][]RankTensorElement = initializeRankTensor(p)
	var currRankVector = make([]int, len(p.Mat[0]))//intialize to all zeros, representing the consensus sequence.
	var consensusKey uint64 = dnaTwoBit.BasesToUint64RightAln(currSeq.Seq, 0, len(currSeq.Seq))
	answer[consensusKey] = consensusValue

	for column := 0; column < len(rankMatrix[0]); column++ {
		//make child sequence and rank vector
		currSeq.Seq[column] = rankMatrix[1][column].Base
		currRankVector[column] = 1
		//call recursive part
		recursiveCheckKmers(answer, currSeq.Seq, rankMatrix, consensusValue, currRankVector, column, threshold)
		//revert child sequence to parent before next recursive call
		currSeq.Seq[column] = rankMatrix[0][column].Base
		currRankVector[column] = 0
	}

	return answer
}

func recursiveCheckKmers(answer map[uint64]float64, currSeq []dna.Base, rankMatrix [][]RankTensorElement, parentValue float64, rankVector []int, index int, threshold float64) {
	//first score sequence
	//we look at parent score, and we look at changed base.
	currValue := parentValue + rankMatrix[rankVector[index]][index].Value - rankMatrix[rankVector[index]-1][index].Value
	// if we pass threshold. Everything else is in this if, so if we fail threshold, we just return
	if currValue >= threshold {
		//if we've passed, we add current sequence to the map
		var currKey uint64 = dnaTwoBit.BasesToUint64RightAln(currSeq, 0, len(currSeq))
		answer[currKey] = currValue
		for i := index; i < len(rankMatrix[0]); i++ {//generate children to the right of the current index.
			if rankVector[i] < 3 { //we only have to generate children if the rank vector is less than 3 at this position (only 4 base identities)
				currSeq[i] = rankMatrix[rankVector[i]][i].Base
				rankVector[i]++
				recursiveCheckKmers(answer, currSeq, rankMatrix, currValue, rankVector, i, threshold)
				//decrement before we do next recursive call
				rankVector[i]--
				currSeq[i] = rankMatrix[rankVector[i]][i].Base//this will decrement because we just decremented rankVector.
			}
		}
	}
}
