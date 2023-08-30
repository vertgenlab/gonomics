package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
)

// incrementWindowEdge takes two aligned fasta sequences and a current index into the alignment.
// the current index can be -1 for when we want to "increment" to the first position in the alignment.
// the sequences are assumed to be all uppercase DNA
// the function returns:
//  1. the next alignment position with a reference base, which is usually alignmentIndex+1, but can be greater due to gaps
//     this value will be equal to len(seqOne) if we were unable to find another reference base
//  2. gaps opened in ref (0 or 1)
//  3. gaps opened in the (query) second sequence (0 or 1)
//  4. gaps closed in the (query) second sequence (0 or 1)
//  5. int that is the number of Ns in the reference after this increment (0 or 1)
//  6. int that is the number of Ns in the query (second) sequence after this increment from matches
//  7. int that is the number of Ns in the query (second) sequence after this increment from gaps in the reference
//  8. int that is the number of substitutions/mismatches (0 or 1)
func incrementWindowEdge(reference []dna.Base, firstQuery []dna.Base, secondQuery []dna.Base, alnIdxOrig int) (alnIdx, gapOpenCloseRef, gapOpenedQuery, gapClosedQuery, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst int) {
	alnIdx = alnIdxOrig

	// increment alnIdx to the next ref position (next non-gap pos in reference)
	for alnIdx++; alnIdx < len(reference) && reference[alnIdx] == dna.Gap; alnIdx++ {
		if seqTwo[alnIdx] == dna.N {
			numQueryNsGap++
		}
		if seqTwo[alnIdx] != dna.Gap {
			gapOpenCloseRef = 1
		}
	}

	// if we ran off the end of seqOne when looking for the next ref base
	if alnIdx == len(seqOne) {
		return
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

	return
}

// speedyWindowDifference is a helper function of faFindFast that calculates the divergence between two input sequences for every position using a sliding window.
// optional arguments longOutput and divergenceRate allow the user to report a -log10pValue corresponding to the p value of observing a level of divergence for a given
// window under a null binomial model of neutral evolution.
func speedyWindowDifference(reference []dna.Base, firstQuery []dna.Base, secondQuery []dna.Base, s Settings) {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1                                                   // these are the two edges of the sliding window in "alignment coordinates"
	var refIdxBeforeWindow, lastRefIdxOfWindow int = -1, -1                                                   // these are the two edges of the sliding window in "reference (no gaps) coordinates"
	var totalGaps, totalNs, totalSubst int                                                                    // this is the data we need to keep track of that describes the current window
	var gapOpenCloseRef, gapOpenQuery, gapClosedQuery, numRefNs, numQueryNsGap, numQueryNsMatch, numSubst int // ints we will get back when moving the window one ref base. TODO: pos ref?
	var err error
	var percentDiverged, rawPValue float64

	// this caches map[k] to -log10(BinomialDist(n, k, p, true)), which is the -log10 p Value.
	var scorePValueCache map[int]float64

	// create outFile
	file := fileio.EasyCreate(s.OutFile)

	// divergenceRate = -1 is a reserved value that signifies that the user has not set a divergence rate. If divergenceRate != -1,
	// we initialize the scorePValueCache.
	if s.DivergenceRate != -1 {
		scorePValueCache = binomialDistCacheLog10(s.WindowSize, s.DivergenceRate)
	}

	for lastAlnIdxOfWindow < len(firstQuery) { // this check could also be "!done", I am not sure what is more clear
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
				if longOutput {
					percentDiverged = 100 * (float64(totalSubst+totalGaps) / float64(windowSize))
					if totalSubst+totalGaps > windowSize {
						log.Fatalf("Error: total number of mutations exceeds windowSize. This may or may not be a bug, but your sequence has deviated from our use case.")
					}
					rawPValue = scorePValueCache[totalSubst+totalGaps]
					_, err = fmt.Fprintf(out, "%s\t%d\t%d\t%s_%d\t%d\t%s\t%e\t%e\n", posRefChrName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, posRefChrName, refIdxBeforeWindow+1, totalSubst+totalGaps, "+", percentDiverged, rawPValue)
					exception.PanicOnErr(err)
				} else {
					_, err = fmt.Fprintf(out, "%s\t%d\t%d\t%s_%d\t%d\n", posRefChrName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, posRefChrName, refIdxBeforeWindow+1, totalSubst+totalGaps)
					exception.PanicOnErr(err)
				}
			}
		}
	}

	// close outFile
	err = file.Close()
	exception.PanicOnErr(err)
}

// binomialDistCacheLog10 generates a map of form map[int]float64 for a binomial distribution specified by parameters n and k.
// map[k] returns a float64 representing -log10(BinomialDist(n, k, p, true)), or the -log10pValue.
func binomialDistCacheLog10(n int, p float64) map[int]float64 {
	if p < 0 || p > 1 {
		log.Fatalf("Error: p must be a value between 0 and 1. Found: %v.\n", p)
	}
	var answer = make(map[int]float64)
	answer[n] = numbers.BinomialDistLog(n, n, p)

	for k := n - 1; k >= 0; k-- {
		answer[k] = logspace.Add(numbers.BinomialDistLog(n, k, p), answer[k+1])
	}

	// now convert to -log10 space
	for k := 0; k <= n; k++ {
		answer[k] = -1 * logspace.ToBase10(answer[k])
	}
	answer[0] = 0 //hardcoded to avoid numerical noise

	return answer
}
