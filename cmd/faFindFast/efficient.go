package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func efficientWindowDifference(windowSize int, seq1 fasta.Fasta, seq2 fasta.Fasta, name *string, verbose bool) []bed.Bed {
	var bedList []bed.Bed
	var alignmentStart, referenceCounter int
	var reachedEnd bool = false
	var containsN bool = false
	var prevNoN bool = false
	var lostGapStart bool = false
	var newGapEnd bool = false
	var gapCount int = 0
	var diff int = 0
	var Ndex int

	for alignmentCounter := 0; reachedEnd == false && (alignmentCounter+windowSize+gapCount) < len(seq1.Seq); alignmentCounter++ {
		if alignmentCounter%1000000 == 0 && verbose {
			fmt.Printf("alignmentCounter: %v\n", alignmentCounter)
		}

		if prevNoN {
			lostGapStart = false
			newGapEnd = false

			//find new referenceStart
			for alignmentStart < len(seq1.Seq) && seq1.Seq[alignmentStart] == dna.Gap {
				alignmentStart++
				gapCount--
				lostGapStart = true
			}
			//find new referenceEnd
			for (alignmentStart+windowSize+gapCount) < len(seq1.Seq) && seq1.Seq[alignmentStart+windowSize+gapCount] == dna.Gap {
				gapCount++
				newGapEnd = true
			}

			if lostGapStart {
				diff--
			}

			if newGapEnd {
				diff++
			}

			if (alignmentStart + windowSize + gapCount) < len(seq1.Seq) {
				if seq1.Seq[alignmentCounter+windowSize] == dna.N || seq2.Seq[alignmentCounter+windowSize] == dna.N {
					prevNoN = false
					alignmentCounter = alignmentCounter + windowSize - 1
					referenceCounter = referenceCounter + windowSize - gapCount - 1
					continue
					//check for added diff at the end
				} else if seq1.Seq[alignmentCounter+windowSize+gapCount] != seq2.Seq[alignmentCounter+windowSize+gapCount] {
					if seq2.Seq[alignmentCounter+windowSize+gapCount] == dna.Gap && seq2.Seq[alignmentCounter+windowSize+gapCount-1] == dna.Gap {
						//do nothing, indel overhang case
					} else {
						diff++
					}
					//check for removed diff at the start
				} else if seq1.Seq[alignmentCounter] != seq2.Seq[alignmentCounter] {
					if seq2.Seq[alignmentCounter] == dna.Gap && seq2.Seq[alignmentCounter-1] == dna.Gap {
						//do nothing, indel overhang case
					} else {
						diff--
					}
				}
				//generate output bed and add to the list
				current := bed.Bed{Chrom: *name, ChromStart: alignmentCounter, ChromEnd: alignmentCounter + windowSize, Name: fmt.Sprintf("%d", referenceCounter), Score: diff, FieldsInitialized: 5}
				bedList = append(bedList, current)
				referenceCounter++
			} else {
				reachedEnd = true
			}
		} else if seq1.Seq[alignmentCounter] != dna.Gap {
			containsN, Ndex, gapCount = efficientSeqsContainN(seq1, seq2, alignmentCounter, windowSize)
			if containsN {
				alignmentCounter = Ndex
				referenceCounter = referenceCounter + Ndex - gapCount //check with Craig
			} else {
				prevNoN = true
				diff, gapCount, reachedEnd = efficientCountWindowDifference(seq1, seq2, alignmentCounter, windowSize)
				if !reachedEnd {
					current := bed.Bed{Chrom: *name, ChromStart: referenceCounter,
						ChromEnd: referenceCounter + windowSize, Name: fmt.Sprintf("%d", referenceCounter), Score: diff, FieldsInitialized: 5}
					bedList = append(bedList, current)
					referenceCounter++
				}
			}
		}
	}
	return bedList
}

func efficientSeqsContainN(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) (bool, int, int) {
	var gapCount int = 0
	var i int
	for i = start; i < windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == dna.Gap {
			gapCount++
		}

		if seq1.Seq[i] == dna.N {
			return true, i, gapCount
		} else if seq2.Seq[i] == dna.N {
			return true, i, gapCount
		}
	}
	return false, i, gapCount
}

func efficientCountWindowDifference(seq1 fasta.Fasta, seq2 fasta.Fasta, start int, windowSize int) (int, int, bool) {
	diff := 0
	baseCount := 0
	var gapCount int = 0
	var seq1Indel bool = false
	var seq2Indel bool = false
	var reachedEnd bool = false
	var i int = 0

	for i = start; baseCount < windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
				baseCount++
			} else {
				gapCount++
			}
		} else if seq1.Seq[i] == dna.Gap {
			gapCount++
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			baseCount++
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			baseCount++
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}

	if baseCount != windowSize {
		reachedEnd = true
	}
	return diff, gapCount, reachedEnd
}
