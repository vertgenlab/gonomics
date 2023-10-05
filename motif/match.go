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
	"math"
)

// ScoreWindow returns the motif score for an input PositionMatrix for a given []dna.Base at a specified start position alnStart.
// First return (float64) is the motif score. Second return (int) is the end alnPos index of the window.
// Third return (bool) returns false
// if the sequence could not be scored (contains Ns, ran off end of sequence).
func ScoreWindow(pm PositionMatrix, seq []dna.Base, alnStart int) (float64, int, bool) {
	var currAlnPos = alnStart
	var answer float64
	var motifPos int = 0
	var seenAnswer string = ""
	for motifPos < len(pm.Mat[0]) {
		if currAlnPos >= len(seq) {
			return -1, -1, false //ran off the end of the sequence.
		}
		switch seq[currAlnPos] {
		case dna.Gap:
			//we advance in alnPos at the end of the loop to check the next position but do not advance the motifPosition.
		case dna.A:
			answer += pm.Mat[0][motifPos]
			motifPos++
			seenAnswer += "A"
		case dna.C:
			answer += pm.Mat[1][motifPos]
			motifPos++
			seenAnswer += "C"
		case dna.G:
			answer += pm.Mat[2][motifPos]
			motifPos++
			seenAnswer += "G"
		case dna.T:
			answer += pm.Mat[3][motifPos]
			motifPos++
			seenAnswer += "T"
		case dna.N:
			return -1, -1, false //this will just ignore this position of the sequence and move on.
		default:
			log.Fatalf("Unrecognized base. Cannot score window.")
			return -1, -1, false
		}
		currAlnPos++
	}
	return answer, currAlnPos, true
}

// RapidMatch performs genome-wide scans for TF motif occurrences from an input genome in fasta format.
// propMatch specifies the motif score threshold for a match, as a proportion of the consensus score.
// outputAsProportion formats the output score reporting as match proportion of consensus score.
func RapidMatch(motifs []PositionMatrix, records []fasta.Fasta, propMatch float64, outFile string, outputAsProportion bool) {
	var err error
	var motifLen int
	var kmerHash map[uint64]float64
	var consensusScore float64
	var couldScoreConsensus bool
	var revCompMotif PositionMatrix
	out := fileio.EasyCreate(outFile)

	for i := range motifs {
		motifLen = len(motifs[i].Mat[0])
		if motifLen > 32 {
			log.Fatalf("RapidMatch cannot accommodate Position Matrices with a motif length greater than 32. Plese filter out the matrix with this ID: %v.\n", motifs[i].Id)
		}
		var currSeq = ConsensusSequence(motifs[i], false)
		consensusScore, _, couldScoreConsensus = ScoreWindow(motifs[i], currSeq.Seq, 0)
		if !couldScoreConsensus {
			log.Fatalf("Error in buildKmerHash. Could not score consensus sequence.")
		}
		kmerHash = buildKmerHash(motifs[i], propMatch)
		scanGenome(records, kmerHash, consensusScore, motifs[i].Name, motifLen, out, bed.Positive, outputAsProportion)
		revCompMotif = ReverseComplement(motifs[i])
		kmerHash = buildKmerHash(revCompMotif, propMatch)
		scanGenome(records, kmerHash, consensusScore, motifs[i].Name, motifLen, out, bed.Negative, outputAsProportion)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

// scanGenome is a helper function of rapid match that scans the genome for motif occurrences for an individual motif,
// whose kmer match scores are embedded in an input hash.
func scanGenome(records []fasta.Fasta, kmerHash map[uint64]float64, consensusScore float64, motifName string, motifLen int, out *fileio.EasyWriter, strand bed.Strand, outputAsProportion bool) {
	var currChrom, currPos, newKeyPos int
	var needNewKey, couldGetNewKey, inKmerHash bool
	var currBed bed.Bed
	var currKey uint64
	var currScore float64
	var bitMask uint64 = uint64(math.Pow(2, float64(2*motifLen)) - 1) // bitmask formula: B_n = 2^{2n} - 1

	for currChrom = range records {
		needNewKey = true // we'll need a new key when we start each chromosome
		for currPos = 0; currPos < len(records[currChrom].Seq); currPos++ {
			if needNewKey {
				currKey, newKeyPos, couldGetNewKey = getNewKey(records[currChrom], currPos, motifLen)
				currPos = newKeyPos
				if !couldGetNewKey {
					break // this means we've run out of windows on the current chrom and we should move to the next chrom.
				}
				needNewKey = false
			} else {
				switch records[currChrom].Seq[currPos] {
				case dna.N:
					needNewKey = true
					continue
				case dna.Gap:
					continue
				case dna.A:
					currKey = currKey << 2
					currKey = currKey | 0
					currKey = currKey & bitMask
				case dna.C:
					currKey = currKey << 2
					currKey = currKey | 1
					currKey = currKey & bitMask
				case dna.G:
					currKey = currKey << 2
					currKey = currKey | 2
					currKey = currKey & bitMask
				case dna.T:
					currKey = currKey << 2
					currKey = currKey | 3
					currKey = currKey & bitMask
				default:
					log.Fatalf("Unrecognized base: %s.", dna.BaseToString(records[currChrom].Seq[currPos]))
				}
			}
			if !needNewKey {
				if currScore, inKmerHash = kmerHash[currKey]; inKmerHash { //if we get a hit from the current key in the kmer hash.
					if outputAsProportion {
						currScore = currScore / consensusScore
					}
					currBed = bed.Bed{Chrom: records[currChrom].Name,
						ChromStart:        currPos - motifLen,
						ChromEnd:          currPos,
						Name:              motifName,
						Score:             0,
						Strand:            strand,
						FieldsInitialized: 7,
						Annotation:        []string{fmt.Sprintf("%f", currScore)}}
					bed.WriteBed(out, currBed)
				}
			}
		}
	}
}

// getNewKey generates a 2bit encoded kmer sequence as a uint64 from a multiple alignment starting at an input alnPos.
// for a pairwise alignment, index should be 0 for target and 1 for query.
// returns the uint64 seq, an int corresponding to the ending alnPos, and a bool that returns true if a key was generated without
// reaching the end of the alignment.
func getNewKey(record fasta.Fasta, alnPos int, motifLen int) (uint64, int, bool) {
	var answer uint64
	var motifPos = 0
	for motifPos < motifLen {
		if alnPos >= len(record.Seq) {
			return 0, 0, false
		}
		switch record.Seq[alnPos] {
		case dna.N: //if we see an N, we skip this base and clear out answer and the partial key.
			motifPos = 0
			answer = 0
		case dna.Gap: //keep searching, motifs can span gaps, which are ignored.
		case dna.A:
			answer = answer << 2 //left shift two bases, clearing a space for the new base.
			answer = answer | 0  //append 00 in last two spots with bitwise OR.
			motifPos++
		case dna.C:
			answer = answer << 2
			answer = answer | 1
			motifPos++
		case dna.G:
			answer = answer << 2
			answer = answer | 2
			motifPos++
		case dna.T:
			answer = answer << 2
			answer = answer | 3
			motifPos++
		default:
			log.Fatalf("Unrecognized base: %s.", dna.BaseToString(record.Seq[alnPos]))
		}
		alnPos++
	}
	return answer, alnPos, true
}

// rankTensorElement encodes the value and base for a rankTensor at a particular row and column.
// In three dimensions, it is useful to consider the Value and Base as two layers of the rankTensor.
type rankTensorElement struct {
	Value float64
	Base  dna.Base
}

// initializeRankTensor turns a PositionMatrix into a rank-ordered position weight tensor.
func initializeRankTensor(p PositionMatrix) [][]rankTensorElement {
	var row, column int
	var currMaxRow, currRank int
	var currMaxValue float64
	var answer = make([][]rankTensorElement, 4)
	for row = 0; row < 4; row++ {
		answer[row] = make([]rankTensorElement, len(p.Mat[row]))
		for column = 0; column < len(p.Mat[row]); column++ {
			answer[row][column] = rankTensorElement{Value: p.Mat[row][column], Base: dna.Base(row)} //cast row to dna.Base (0 -> A, 1 -> C, 2 -> G, 3 -> T)
		}
	}
	//sort matrix columns by rank
	for column = 0; column < len(p.Mat[0]); column++ {
		for currRank = 0; currRank < 3; currRank++ {
			currMaxRow = currRank
			currMaxValue = answer[currRank][column].Value
			for row = currRank + 1; row < 4; row++ {
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

// rankTensorToString formats a rankTensor as a string for debugging and visualization.
func rankTensorToString(m [][]rankTensorElement) string {
	var row, column int
	var answer, currBaseString string
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

// buildKmerHash produces a hash mapping 2bit encoded kmer sequences to their corresponding motif score for a PositionMatrix.
// Only kmers with a motif score above an input thresholdProportion are stored in the output map.
func buildKmerHash(p PositionMatrix, thresholdProportion float64) map[uint64]float64 {
	var answer = make(map[uint64]float64)
	var currSeq = ConsensusSequence(p, false)
	consensusValue, _, couldScoreConsensus := ScoreWindow(p, currSeq.Seq, 0)
	if !couldScoreConsensus {
		log.Fatalf("Error in buildKmerHash. Could not score consensus sequence.")
	}
	var threshold float64 = thresholdProportion * consensusValue
	var rankMatrix [][]rankTensorElement = initializeRankTensor(p)
	var currRankVector = make([]int, len(p.Mat[0])) //initialize to all zeros, representing the consensus sequence.
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

// recursiveCheckKmers is a helper function of buildKmerHash, which searches for all kmer permutations of a PostionMatrix consensus sequence
// and adds all kmers with a motif score above an input threshold to a map. This implementation uses dynamic programming to search all kmers efficiently.
func recursiveCheckKmers(answer map[uint64]float64, currSeq []dna.Base, rankMatrix [][]rankTensorElement, parentValue float64, rankVector []int, index int, threshold float64) {
	//first score sequence
	//we look at parent score, and we look at changed base.
	currValue := parentValue + rankMatrix[rankVector[index]][index].Value - rankMatrix[rankVector[index]-1][index].Value
	// if we pass threshold. Everything else is in this if, so if we fail threshold, we just return
	if currValue >= threshold {
		//if we've passed, we add current sequence to the map
		var currKey uint64 = dnaTwoBit.BasesToUint64RightAln(currSeq, 0, len(currSeq))
		answer[currKey] = currValue
		for i := index; i < len(rankMatrix[0]); i++ { //generate children to the right of the current index.
			if rankVector[i] < 3 { //we only have to generate children if the rank vector is less than 3 at this position (only 4 base identities)
				currSeq[i] = rankMatrix[rankVector[i]][i].Base
				rankVector[i]++
				recursiveCheckKmers(answer, currSeq, rankMatrix, currValue, rankVector, i, threshold)
				//decrement before we do next recursive call
				rankVector[i]--
				currSeq[i] = rankMatrix[rankVector[i]][i].Base //this will decrement because we just decremented rankVector.
			}
		}
	}
}
