package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"io"
)

// incrementWindowEdge takes two aligned fasta sequences and a current index into the alignment.
// the current index can be -1 for when we want to "increment" to the first position in the alignment.
// the sequences are assumed to be all uppercase DNA
// the function returns:
// 1) the next alignment position with a reference base, which is usually alignmentIndex+1, but can be greater due to gaps
//        this value will be equal to len(seqOne) if we were unable to find another reference base
// 2) gaps opened in ref (0 or 1)
// 3) gaps opened in the (query) second sequence (0 or 1)
// 4) gaps closed in the (query) second sequence (0 or 1)
// 5) int that is the number of Ns in the reference after this increment (0 or 1)
// 6) int that is the number of Ns in the query (second) sequence after this increment from matches
// 7) int that is the number of Ns in the query (second) sequence after this increment from gaps in the reference
// 8) int that is the number of substitutions/mismatches (0 or 1)
func incrementWindowEdge(seqOne []dna.Base, seqTwo []dna.Base, alnIdx int) (int, int, int, int, int, int, int, int) {
	var alnIdxOrig, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst, gapOpenCloseRef, gapOpenedQuery, gapClosedQuery int
	alnIdxOrig = alnIdx

	// increment alnIdx to the next ref position (next non-gap pos in seqOne)
	for alnIdx++; alnIdx < len(seqOne) && seqOne[alnIdx] == dna.Gap; alnIdx++ {
		if seqTwo[alnIdx] == dna.N {
			numQueryNsGap++
		}
		if seqTwo[alnIdx] != dna.Gap {
			gapOpenCloseRef = 1
		}
	}

	// if we ran off the end of seqOne when looking for the next ref base
	if alnIdx == len(seqOne) {
		return alnIdx, 0, 0, 0, 0, numQueryNsGap, 0, 0
	}

	// did we add another reference N when moving the window one reference base
	if seqOne[alnIdx] == dna.N {
		numRefNs++
	}
	// do we add another N to the query count of Ns
	if seqTwo[alnIdx] == dna.N {
		numQueryNsMatch++
	}
	// is this a substitution?
	if seqOne[alnIdx] != seqTwo[alnIdx] && dna.DefineBase(seqOne[alnIdx]) && dna.DefineBase(seqTwo[alnIdx]) {
		numSubst++
	}
	// did we open a gap in the query sequence when moving the window edge?
	if ((alnIdxOrig != -1 && seqTwo[alnIdxOrig] != dna.Gap) || alnIdxOrig == -1) && seqTwo[alnIdx] == dna.Gap {
		gapOpenedQuery++
	}
	// did we close a gap in the query when moving the window edge?
	if alnIdxOrig != -1 && seqTwo[alnIdxOrig] == dna.Gap && seqTwo[alnIdx] != dna.Gap {
		//if seqTwo[alnIdx] == dna.Gap && (alnIdx+1==len(seqOne) || seqTwo[alnIdx+1] != dna.Gap) {
		gapClosedQuery++
	}

	return alnIdx, gapOpenCloseRef, gapOpenedQuery, gapClosedQuery, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst
}

func speedyWindowDifference(windowSize int, reference []dna.Base, query []dna.Base, refChrName string, noPrintIfN bool, verbose bool, out io.Writer) {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1                                                   // these are the two edges of the sliding window in "alignment coordinates"
	var refIdxBeforeWindow, lastRefIdxOfWindow int = -1, -1                                                   // these are the two edges of the sliding window in "reference (no gaps) coordinates"
	var totalGaps, totalNs, totalSubst int                                                                    // this is the data we need to keep track of that describes the current window
	var gapOpenCloseRef, gapOpenQuery, gapClosedQuery, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst int // ints we will get back when moving the window one ref base

	for lastAlnIdxOfWindow < len(reference) { // this check could also be "!done", I am not sure what is more clear
		// we always move the lastBaseOfTheWindow (right side) and add on what we find to the counters because
		// all this stuff is now inside the current window
		lastAlnIdxOfWindow, gapOpenCloseRef, gapOpenQuery, _, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst = incrementWindowEdge(reference, query, lastAlnIdxOfWindow)
		lastRefIdxOfWindow++
		totalGaps += gapOpenCloseRef + gapOpenQuery
		totalNs += numRefNs + numQueryNsGap + numQueryNsMatch
		totalSubst += numSubst

		// usually increment the baseBeforeWindow, but not at the beginning when we have not yet incremented the end
		// enough to have a full "windowSize" of bases in the window
		if lastRefIdxOfWindow-refIdxBeforeWindow > windowSize {
			alnIdxBeforeWindow, _, _, _, numRefNs, _, numQueryNsMatch, numSubst = incrementWindowEdge(reference, query, alnIdxBeforeWindow)
			refIdxBeforeWindow++
			//totalGaps -= gapOpenCloseRef + gapClosedQuery
			totalNs -= numRefNs + numQueryNsMatch
			totalSubst -= numSubst
		}

		// the trailing window needs to "look ahead" to see what happens before the next ref base
		if lastRefIdxOfWindow-refIdxBeforeWindow == windowSize {
			_, gapOpenCloseRef, _, gapClosedQuery, _, numQueryNsGap, _, _ = incrementWindowEdge(reference, query, alnIdxBeforeWindow)
			totalGaps -= gapOpenCloseRef + gapClosedQuery
			totalNs -= numQueryNsGap
		}

		// we check to make sure we are not at the very beginning or end, where we would have partial or illegal windows
		if lastRefIdxOfWindow-refIdxBeforeWindow == windowSize && lastAlnIdxOfWindow < len(reference) {
			// an option/flag can tell us not to print if there are Ns in the query or target
			if !noPrintIfN || totalNs == 0 {
				fmt.Fprintf(out, "%s\t%d\t%d\t%s_%d\t%d\n", refChrName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, refChrName, refIdxBeforeWindow+1, totalSubst+totalGaps)
			}
		}
	}
}
