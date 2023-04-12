package simulate

import (
	"fmt"
	"log"
	"math"
	"math/rand"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
)

//RandIntergenicSeq makes a randomly generated DNA sequence of a specified length and GC content. Unlike RandGene, it does not have to be divisible by 3.
func RandIntergenicSeq(GcContent float64, lenSeq int) []dna.Base {
	var answer []dna.Base = make([]dna.Base, lenSeq)
	for i := range answer {
		answer[i] = chooseRandomBase(GcContent)
	}
	return answer
}

func indelLength(lambda float64) int {
	expFloat, _ := numbers.RandExp()
	return int(math.Ceil(expFloat / lambda))
}

const bufferSize = 10_000_000

// SimulateWithIndels takes an input fastaFile, which must contain to a single fasta entry, and simulates a mutated sequence.
// The output sequence is provided in a multiFa alignment, aligned ot the initial sequence.
// branchLength (a float from 0 to 1) specifies the expected value of the proportion of sites in the input sequence that will be mutated.
// propIndel (a float from 0 to 1)specifies the expected value of the proportion of indels in the output sequence.
// lambda specifies the rate parameter for an exponential distribution, from which simulated INDEL sizes will be sampled.
// gcContent specifies the expected value of GC content for inserted sequences.
// vcfOutFile specifies an optional return, which records all variants made during the simulated mutation process.
// transitionBias specifies the expected value of the ratio of transitions to transversions in the output sequence.
// qName sets the suffix for the output query fasta name.
func SimulateWithIndels(fastaFile string, branchLength float64, propIndel float64, lambda float64, gcContent float64, transitionBias float64, vcfOutFile string, qName string) []fasta.Fasta {
	var answer = make([]fasta.Fasta, 2)
	var emptyRoomInBuffer = bufferSize
	var currRand, currRand2, currRand3 float64
	var outputPos, length, indelPos, indelStartPos int
	var newBufferRoom []dna.Base = make([]dna.Base, bufferSize)
	var inputPos int = 0
	var vcfOut *fileio.EasyWriter
	var currRef, currAlt []dna.Base
	var err error
	var outOfChrom bool = false
	records := fasta.Read(fastaFile)
	if len(records) != 1 {
		log.Fatalf("SimulateWithIndels expects a single fasta record in the input file.")
	}
	answer[0] = fasta.Fasta{Name: records[0].Name, Seq: make([]dna.Base, bufferSize)}
	answer[1] = fasta.Fasta{Name: fmt.Sprintf("%s_%s", records[0].Name, qName), Seq: make([]dna.Base, bufferSize)}

	if vcfOutFile != "" {
		vcfOut = fileio.EasyCreate(vcfOutFile)
		_, err = fmt.Fprintf(vcfOut, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
		exception.PanicOnErr(err)
	}

	for inputPos < len(records[0].Seq) {
		currRand = rand.Float64() //this rand determines if there will be a mutation
		if currRand < branchLength {
			currRand2 = rand.Float64()         //this rand determines the mutation type
			if currRand2 < (propIndel / 2.0) { //this half of INDELs will be deletions
				indelStartPos = inputPos + 1
				currRand3 = rand.Float64() //one more rand for the case where a substitution immediately precedes a deletion
				if currRand3 < branchLength {
					answer[0].Seq[outputPos] = records[0].Seq[inputPos]
					currRef = []dna.Base{records[0].Seq[inputPos]}
					if transitionBias != 1 {
						answer[1].Seq[outputPos] = changeBaseTransitionBias(records[0].Seq[inputPos], transitionBias)
					} else {
						answer[1].Seq[outputPos] = changeBase(records[0].Seq[inputPos])
					}
					currAlt = []dna.Base{answer[1].Seq[outputPos]}
				} else {
					answer[0].Seq[outputPos] = records[0].Seq[inputPos]
					currRef = []dna.Base{records[0].Seq[inputPos]}
					answer[1].Seq[outputPos] = records[0].Seq[inputPos]
					currAlt = []dna.Base{records[0].Seq[inputPos]}
				}
				inputPos++
				if inputPos >= len(records[0].Seq) {
					break
				}
				outputPos++
				emptyRoomInBuffer--
				if emptyRoomInBuffer < 1 {
					answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
					answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
					emptyRoomInBuffer += bufferSize
				}
				length = indelLength(lambda)
				indelPos = 0
				for indelPos < length {
					answer[0].Seq[outputPos] = records[0].Seq[inputPos]
					currRef = append(currRef, records[0].Seq[inputPos])
					answer[1].Seq[outputPos] = dna.Gap
					outputPos++
					emptyRoomInBuffer--
					if emptyRoomInBuffer < 1 {
						answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
						answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
						emptyRoomInBuffer += bufferSize
					}
					indelPos++
					inputPos++
					if inputPos >= len(records[0].Seq) {
						outOfChrom = true
						break
					}
				}
				inputPos--      //avoids double skipping after indel
				if outOfChrom { //double break needed if we run out of chromosome in the inner loop
					break
				}
				//if we didn't run off the chrom, we'll report the variant
				if vcfOutFile != "" {
					_, err = fmt.Fprintf(vcfOut, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n", records[0].Name, indelStartPos, ".", dna.BasesToString(currRef), dna.BasesToString((currAlt)), "100", "PASS", ".", ".")
				}
			} else if currRand2 < propIndel { //the other half will be insertions
				indelStartPos = inputPos + 1
				currRand2 = rand.Float64()
				if currRand2 < branchLength { //case where a substitution immediately precedes an insertion
					answer[0].Seq[outputPos] = records[0].Seq[inputPos]
					currRef = []dna.Base{records[0].Seq[inputPos]}
					if transitionBias != 1 {
						answer[1].Seq[outputPos] = changeBaseTransitionBias(records[0].Seq[inputPos], transitionBias)
					} else {
						answer[1].Seq[outputPos] = changeBase(records[0].Seq[inputPos])
					}
					currAlt = []dna.Base{answer[1].Seq[outputPos]}
				} else {
					answer[0].Seq[outputPos] = records[0].Seq[inputPos]
					currRef = []dna.Base{records[0].Seq[inputPos]}
					answer[1].Seq[outputPos] = records[0].Seq[inputPos]
					currAlt = []dna.Base{records[0].Seq[inputPos]}
				}
				inputPos++
				if inputPos >= len(records[0].Seq) {
					break
				}
				outputPos++
				emptyRoomInBuffer--
				if emptyRoomInBuffer < 1 {
					answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
					answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
					emptyRoomInBuffer += bufferSize
				}
				length = indelLength(lambda)
				indelPos = 0
				for indelPos < length {
					answer[0].Seq[outputPos] = dna.Gap
					answer[1].Seq[outputPos] = chooseRandomBase(gcContent)
					currAlt = append(currAlt, answer[1].Seq[outputPos])
					outputPos++
					emptyRoomInBuffer--
					if emptyRoomInBuffer < 1 {
						answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
						answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
						emptyRoomInBuffer += bufferSize
					}
					indelPos++
				}
				inputPos-- //avoid double skipping for the main loop iterator
				if vcfOutFile != "" {
					_, err = fmt.Fprintf(vcfOut, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n", records[0].Name, indelStartPos, ".", dna.BasesToString(currRef), dna.BasesToString(currAlt), "100", "PASS", ".", ".")
				}
			} else { //this section handles the substitution case
				answer[0].Seq[outputPos] = records[0].Seq[inputPos]
				if transitionBias != 1 {
					answer[1].Seq[outputPos] = changeBaseTransitionBias(records[0].Seq[inputPos], transitionBias)
				} else {
					answer[1].Seq[outputPos] = changeBase(records[0].Seq[inputPos])
				}
				currRef = []dna.Base{records[0].Seq[inputPos]}
				currAlt = []dna.Base{answer[1].Seq[outputPos]}
				if vcfOutFile != "" {
					_, err = fmt.Fprintf(vcfOut, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n", records[0].Name, inputPos+1, ".", dna.BasesToString(currRef), dna.BasesToString(currAlt), "100", "PASS", ".", ".")
				}
				outputPos++
				emptyRoomInBuffer--
				if emptyRoomInBuffer < 1 {
					answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
					answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
					emptyRoomInBuffer += bufferSize
				}
			}
		} else {
			answer[0].Seq[outputPos] = records[0].Seq[inputPos]
			answer[1].Seq[outputPos] = records[0].Seq[inputPos]
			outputPos++
			emptyRoomInBuffer--
			if emptyRoomInBuffer < 1 {
				answer[0].Seq = append(answer[0].Seq, newBufferRoom...)
				answer[1].Seq = append(answer[1].Seq, newBufferRoom...)
				emptyRoomInBuffer += bufferSize
			}
		}
		inputPos++
	}

	if vcfOutFile != "" {
		err = vcfOut.Close()
		exception.PanicOnErr(err)
	}

	//trim extra room in sequence from the buffer
	answer[0].Seq = answer[0].Seq[:len(answer[0].Seq)-emptyRoomInBuffer]
	answer[1].Seq = answer[1].Seq[:len(answer[1].Seq)-emptyRoomInBuffer]

	return answer
}

// changeBaseTransitionBias substitutes an input base b following the K80 model with a transitionBias parameter gamma.
func changeBaseTransitionBias(b dna.Base, gamma float64) dna.Base {
	var rand = rand.Float64()
	var TvProb float64 = 1.0 / (2.0 + gamma) //p(a -> c) = 1 / (2 + gamma) (transversion)
	switch dna.ToUpper(b) {
	case dna.A:
		if rand < TvProb { //p(a -> c) = 1 / (2 + gamma) (transversion)
			return dna.C
		} else if rand < 2.0*TvProb { //p(a -> t) = 1 / (2+gamma) (transversion)
			return dna.T
		} else { //that leaves p(a -> g) = gamma / (2+gamma) of probability left, (transition) and p(a -> g) = gamma * p(a -> c) = gamma* p(a -> t)
			return dna.G
		}
	case dna.C:
		if rand < 1.0/(2.0+gamma) {
			return dna.A //transversion
		} else if rand < 2.0/(2.0+gamma) {
			return dna.G //transversion
		} else {
			return dna.T //transition
		}
	case dna.G:
		if rand < 1.0/(2.0+gamma) {
			return dna.C //transversion
		} else if rand < 2.0/(2.0+gamma) {
			return dna.T //transversion
		} else {
			return dna.A //transition
		}
	case dna.T:
		if rand < 1.0/(2.0+gamma) {
			return dna.A //transversion
		} else if rand < 2.0/(2.0+gamma) {
			return dna.G //transversion
		} else {
			return dna.C //transition
		}
	case dna.N:
		return dna.N
	}
	log.Fatalf("Unrecognized base: %v.\n", dna.BaseToString(b))
	return dna.N
}
