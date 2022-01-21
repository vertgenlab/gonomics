package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// incrementWindowEdge takes two aligned fasta sequences and a current index into the alignment.
// the current index can be -1 for when we want to "increment" to the first position in the alignment.
// the sequences are assumed to be all uppercase DNA
// the function returns:
// 1) the next alignment position with a reference base, which is usually alignmentIndex+1, but can be greater due to gaps
//        this value will be equal to len(seqOne) if we were unable to find another reference base
// 2) true if we opened and closed a gap in the ref
// 3) true if a gap was opened in the (query) second sequence
// 4) true if gap was closed in the (query) second sequence
// 5) int that is the number of Ns in the reference after this increment (0 or 1)
// 6) int that is the number of Ns in the query (second) sequence after this increment (can be more than 1 due to gaps in ref)
// 7) int that is the number of substitutions/mismatches (0 or 1)
func incrementWindowEdge(seqOne []dna.Base, seqTwo []dna.Base, alnIdx int) (int, bool, bool, bool, int, int, int) {
	var alnIdxOrig, nextRefPos, numRefNs, numQueryNs, numSubst int
	var gapOpenCloseRef, gapOpenedQuery, gapClosedQuery bool
	alnIdxOrig = alnIdx

	// increment alnIdx to the next ref position (next non-gap pos in seqOne)
	for alnIdx++; alnIdx < len(seqOne) && seqOne[alnIdx] == dna.Gap; alnIdx++ {
		if seqTwo[alnIdx] == dna.Gap {
			numQueryNs++
		}
		if seqTwo[alnIdx] != dna.Gap {
			gapOpenCloseRef = true
		}
	}

	// if we ran off the end of seqOne when looking for the next ref base
	if alnIdx == len(seqOne) {
		return alnIdx, false, false, false, 0, numQueryNs, 0
	}

	// did we add another reference N when moving the window one reference base
	if seqOne[alnIdx] == dna.N {
		numRefNs++
	}
	// do we add another N to the query count of Ns
	if seqTwo[alnIdx] == dna.N {
		numQueryNs++
	}
	// is this a substitution?
	if seqOne[alnIdx] != seqTwo[alnIdx] && dna.DefineBase(seqOne[alnIdx]) && dna.DefineBase(seqTwo[alnIdx]) {
		numSubst++
	}
	// did we open a gap in the query sequence when moving the window edge?
	if alnIdxOrig != -1 && seqTwo[alnIdxOrig] != dna.Gap && seqTwo[alnIdx] == dna.Gap {
		gapOpenedQuery = true
	}
	// did we close a gap in the query when moving the window edge?
	if alnIdxOrig != -1 && seqTwo[alnIdxOrig] == dna.Gap && seqTwo[alnIdx] != dna.Gap {
		gapClosedQuery = true
	}
	return alnIdx, gapOpenCloseRef, gapOpenQuery, gapClosedQuery, numRefNs, numQueryNs, numSubst
}

func speedyWindowDifference(windowSize int, seqOne fasta.Fasta, seqTwo fasta.Fasta, name *string, verbose bool) []bed.Bed {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1 // these are the two edges of the sliding window in "alignment coordinates"
	var refIdxBeforeWindow, lastRefIdxOfWindow int = -1, -1 // these are the two edges of the sliding window in "reference (no gaps) coordinates"
	var totalGaps, totalNs, totalSubst int // this is the data we need to keep track of that describes the current window
	var gapOpenCloseRef, gapOpenQuery, isSubst, gapClosedQuery, done bool // bools we will get back when moving the window one ref base
	var numRefNs, numQueryNs int // ints we will get back when moving the window one ref base

	for lastAlnIdxOfWindow < len(seqOne.Seq) { // this check could also be "!done", I am not sure what is more clear
		// we always move the lastBaseOfTheWindow (right side) and add on what we find to the counters because
		// all this stuff is now inside the current window
		lastAlnIdxOfWindow, gapOpenCloseRef, gapOpenQuery, _, numRefNs, numQueryNs, numSubst = incrementWindowEdge(seqOne.Seq, seqTwo.Seq, lastAlnIdxOfWindow)
		lastRefIdxOfWindow++
		totalGaps += gapOpenCloseRef + gapOpenQuery
		totalNs += numRefNs + numQueryNs
		totalSubst += numSubst

		// usually increment the baseBeforeWindow, but not at the beginning when we have not yet incremented the end
		// enough to have a full "windowSize" of bases in the window
		if lastRefIdxOfWindow - refIdxBeforeWindow < windowSize {
			alnIdxBeforeWindow, gapOpenCloseRef, _, gapClosedQuery, numRefNs, numQueryNs, numSubst = incrementWindowEdge(seqOne.Seq, seqTwo.Seq, alnIdxBeforeWindow)
			refIdxBeforeWindow++
			totalGaps -= gapOpenCloseRef + gapClosedQuery
			totalNs -= numRefNs + numQueryNs
			totalSubst -= numSubst
		}

		// we check to make sure we are not at the very beginning or end, where we would have partial or illegal windows
		if lastRefIdxOfWindow - refIdxBeforeWindow == windowSize && lastAlnIdx < len(seqOne.Seq) {
			fmt.Fprintf()
		}
	}
}

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
