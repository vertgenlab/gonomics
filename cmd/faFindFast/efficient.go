package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
)

// incrementWindowEdge
// inputs: 2 aligned fasta sequences (firstQuery, secondQuery), a current alignment position index.
// the sequences are assumed to be all uppercase DNA
// firstQuery and secondQuery are of equal length
// note that alignment positions keep track of firstQuery and secondQuery sequences, while reference positions keep track of the reference sequence
// outputs:
//  1. alnIdx: the next alignment position where a firstQuery base exists, which is usually alignmentIndex+1, but can be greater due to gaps in the firstQuery sequence
//     this value will be equal to len(firstQuery) if we were unable to find another firstQuery base, and this would mean we have reached the end of the firstQuery sequence
//  2. gapOpenCloseFirstQuery: gaps opened in the firstQuery sequence (0 or 1)
//  3. gapOpenedSecondQuery: gaps opened in the secondQuery sequence (0 or 1)
//  4. gapClosedSecondQuery: gaps closed in the secondQuery sequence (0 or 1)
//  5. numFirstQueryNs: the number of Ns in the firstQuery after this increment (0 or 1)
//  6. numSecondQueryNsGap: the number of Ns in the secondQuery sequence after this increment from gaps in the firstQuery
//  7. numSecondQueryNsMatch: the number of Ns in the secondQuery sequence after this increment from matches with the firstQuery
//  8. numSubst: the number of substitutions/mismatches (0 or 1) between firstQuery and secondQuery
func incrementWindowEdge(firstQuery []dna.Base, secondQuery []dna.Base, alnIdxOrig int) (alnIdx, gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst int) {
	alnIdx = alnIdxOrig
	// at initialization, lastAlnIdxOfWindow == -1, which is assigned to alnIdxOrig when incrementWindowEdge is called, so alnIdx = alnIdxOrig == -1. For loop will start at alnIdx++ which is 0, the 1st position in the alignment

	// increment alnIdx to the next non-gap position in the firstQuery sequence
	// for loop start condition: start at alnIdx++, which is the position after the current alignment position
	// for loop end condition: ends when we reach the end of the firstQuery sequence, or when we reach a non-gap position in the firstQuery sequence, whichever happens first
	// the for loop's start and end conditions also mean that the for loop is never executed, unless there is firstQuery gap before we reach the end of firstQuery
	// for loop incrementation: alnIdx++
	for alnIdx++; alnIdx < len(firstQuery) && firstQuery[alnIdx] == dna.Gap; alnIdx++ {
		// keep track of changes between the original and next non-gap positions in the firstQuery sequence
		// 6. secondQuery
		if secondQuery[alnIdx] == dna.N {
			numSecondQueryNsGap++
		}
		// 2. gapOpenClosedFirstQuery
		if secondQuery[alnIdx] != dna.Gap {
			gapOpenCloseFirstQuery = 1 // Rationale: if the for loop is executed, this means there is a gap in the firstQuery, and if secondQuery does not have a gap, then it is indeed a firstQuery gap compared to secondQuery
		}
	}

	// if we ran off the end of firstQuery when looking for the next non-gap firstQuery base, aka the for loop ended because we reached the end of the firstQuery sequence
	if alnIdx == len(firstQuery) {
		return // return all 8 named return variables
	}

	// 5. numFirstQueryNs: did we add another firstQuery N after moving the window edge by one firstQuery base
	if firstQuery[alnIdx] == dna.N {
		numFirstQueryNs++
	}
	// 7. numSecondQueryNsMatch: did we add another N to the secondQuery count of Ns after moving the window edge by one firstQuery base
	if secondQuery[alnIdx] == dna.N {
		numSecondQueryNsMatch++ // Rationale: these secondQueryNs align with non-gap firstQuery positions
	}
	// 8. numSubst: is this a substitution?
	if firstQuery[alnIdx] != secondQuery[alnIdx] && dna.DefineBase(firstQuery[alnIdx]) && dna.DefineBase(secondQuery[alnIdx]) {
		numSubst++
	}
	// 3. gapOpenedSecondQuery: did we open a gap in the secondQuery sequence when moving the window edge?
	if ((alnIdxOrig != -1 && secondQuery[alnIdxOrig] != dna.Gap) || alnIdxOrig == -1) && secondQuery[alnIdx] == dna.Gap {
		gapOpenedSecondQuery++
	}
	// 4. gapClosedSecondQuery: id we close a gap in the secondQuery when moving the window edge?
	if alnIdxOrig != -1 && secondQuery[alnIdxOrig] == dna.Gap && secondQuery[alnIdx] != dna.Gap {
		//if seqTwo[alnIdx] == dna.Gap && (alnIdx+1==len(seqOne) || seqTwo[alnIdx+1] != dna.Gap) {
		gapClosedSecondQuery++
	}

	return
}

// speedyWindowDifference is a helper function of faFindFast that calculates the divergence between two input sequences (firstQuery, secondQuery sequences) using a sliding window, and then reports the divergence in terms of reference positions (positions in the reference sequence).
// optional arguments longOutput and divergenceRate allow the user to report a -log10pValue corresponding to the p value of observing a level of divergence for a given
// window under a null binomial model of neutral evolution.
func speedyWindowDifference(reference []dna.Base, firstQuery []dna.Base, secondQuery []dna.Base, s Settings) {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1                                                                                           // these are the two edges of the sliding window in "alignment positions"
	var firstQueryIdxBeforeWindow, lastFirstQueryIdxOfWindow int = -1, -1                                                                             // these are the two edges of the sliding window in "firstQuery (no gaps) positions"
	var refIdxBeforeWindow, lastRefIdxOfWindow int = -1, -1                                                                                           // these are the two edges of the sliding window in "reference (no gaps) positions"
	var alnIdxForRef int = 0                                                                                                                          // TODO: remove this variable if needed
	var totalGaps, totalNs, totalSubst int                                                                                                            // this is the data we need to keep track of that describes the current window
	var gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst int // ints we will get back when moving the window one ref base.
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
		lastAlnIdxOfWindow, gapOpenCloseFirstQuery, gapOpenedSecondQuery, _, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst = incrementWindowEdge(firstQuery, secondQuery, lastAlnIdxOfWindow)
		lastFirstQueryIdxOfWindow++
		totalGaps += gapOpenCloseFirstQuery + gapOpenedSecondQuery
		totalNs += numFirstQueryNs + numSecondQueryNsGap + numSecondQueryNsMatch
		totalSubst += numSubst

		// usually increment the baseBeforeWindow,
		// but not at the beginning when we have not yet incremented the end enough to have a full "windowSize" of bases in the window
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow > s.WindowSize {
			alnIdxBeforeWindow, _, _, _, numFirstQueryNs, _, numSecondQueryNsMatch, numSubst = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow)
			firstQueryIdxBeforeWindow++
			//totalGaps -= gapOpenCloseRef + gapClosedQuery
			totalNs -= numFirstQueryNs + numSecondQueryNsMatch
			totalSubst -= numSubst
		}

		// the trailing window needs to "look ahead" to see what happens before the next firstQuery base
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow == s.WindowSize {
			_, gapOpenCloseFirstQuery, _, gapClosedSecondQuery, _, numSecondQueryNsGap, _, _ = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow)
			totalGaps -= gapOpenCloseFirstQuery + gapClosedSecondQuery
			totalNs -= numSecondQueryNsGap
		}

		// we check to make sure we are not at the very beginning or end, where we would have partial or illegal windows
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow == s.WindowSize && lastAlnIdxOfWindow < len(firstQuery) {

			// TODO: convert position from firstQuery to reference, V2 calculate after firstQuery window is finalized, remove all fmt.Printf after debugging

			// V2.1: I calculated refIdx with no other function
			//fmt.Printf("alnIdxBeforeWindow: %v, lastAlnIdxOfWindow: %v, firstQueryIdxBeforeWindow: %v, lastFirstQueryIdxOfWindow: %v\n", alnIdxBeforeWindow, lastAlnIdxOfWindow, firstQueryIdxBeforeWindow, lastFirstQueryIdxOfWindow)

			// check if firstQuery start (non-gap) is gap in reference. If so, don't report that window
			// check if firstQuery end-1 (last non-gap position) is gap in reference. If so, don't report that window
			if !((reference[alnIdxBeforeWindow+1] == dna.Gap && firstQuery[alnIdxBeforeWindow+1] != dna.Gap) || (reference[lastAlnIdxOfWindow] == dna.Gap && firstQuery[lastAlnIdxOfWindow] != dna.Gap)) {

				refIdxBeforeWindow = firstQueryIdxBeforeWindow
				for alnIdxForRef = 0; alnIdxForRef <= alnIdxBeforeWindow; alnIdxForRef++ {
					if reference[alnIdxForRef] == dna.Gap && firstQuery[alnIdxForRef] != dna.Gap {
						//fmt.Printf("alnIdxForRef: %v, reference[alnIdxForRef]: %v, firstQuery[alnIdxForRef]: %v, refIdxBeforeWindow: %v\n", alnIdxForRef, reference[alnIdxForRef], firstQuery[alnIdxForRef], refIdxBeforeWindow)
						refIdxBeforeWindow--
						//fmt.Printf("refIdxBeforeWindow: %v\n", refIdxBeforeWindow)
					} else if reference[alnIdxForRef] != dna.Gap && firstQuery[alnIdxForRef] == dna.Gap {
						//fmt.Printf("alnIdxForRef: %v, reference[alnIdxForRef]: %v, firstQuery[alnIdxForRef]: %v, refIdxBeforeWindow: %v\n", alnIdxForRef, reference[alnIdxForRef], firstQuery[alnIdxForRef], refIdxBeforeWindow)
						refIdxBeforeWindow++
						//fmt.Printf("refIdxBeforeWindow: %v\n", refIdxBeforeWindow)
					}
				}

				lastRefIdxOfWindow = lastFirstQueryIdxOfWindow
				for alnIdxForRef = 0; alnIdxForRef <= lastAlnIdxOfWindow; alnIdxForRef++ {
					if reference[alnIdxForRef] == dna.Gap && firstQuery[alnIdxForRef] != dna.Gap {
						//fmt.Printf("alnIdxForRef: %v, reference[alnIdxForRef]: %v, firstQuery[alnIdxForRef]: %v, lastRefIdxOfWindow: %v\n", alnIdxForRef, reference[alnIdxForRef], firstQuery[alnIdxForRef], lastRefIdxOfWindow)
						lastRefIdxOfWindow--
						//fmt.Printf("lastRefIdxOfWindow: %v\n", lastRefIdxOfWindow)
					} else if reference[alnIdxForRef] != dna.Gap && firstQuery[alnIdxForRef] == dna.Gap {
						//fmt.Printf("alnIdxForRef: %v, reference[alnIdxForRef]: %v, firstQuery[alnIdxForRef]: %v, lastRefIdxOfWindow: %v\n", alnIdxForRef, reference[alnIdxForRef], firstQuery[alnIdxForRef], lastRefIdxOfWindow)
						lastRefIdxOfWindow++
						//fmt.Printf("lastRefIdxOfWindow: %v\n", lastRefIdxOfWindow)
					}
				}

				//fmt.Printf("Will report this window. refIdxBeforeWindow+1:%v, lastRefIdxOfWindow+1:%v\n", refIdxBeforeWindow+1, lastRefIdxOfWindow+1)

				// V2.2: I calculated refIdx with another function fasta.AlnPosToRefPosCounterSeq
				refIdxBeforeWindow_V2 := fasta.AlnPosToRefPosCounterSeq(reference, alnIdxBeforeWindow, 0, 0)
				lastRefIdxOfWindow_V2 := fasta.AlnPosToRefPosCounterSeq(reference, lastAlnIdxOfWindow, 0, 0)
				if refIdxBeforeWindow_V2 != 0 && refIdxBeforeWindow != refIdxBeforeWindow_V2 {
					fmt.Printf("BeforeWindow V1 and V2 don't agree. alnIdxBeforeWindow: %v, refIdxBeforeWindow:%v, refIdxBeforeWindow_V2:%v\n", alnIdxBeforeWindow, refIdxBeforeWindow, refIdxBeforeWindow_V2)
					fmt.Printf("Will report this window. refIdxBeforeWindow+1:%v, lastRefIdxOfWindow+1:%v\n", refIdxBeforeWindow+1, lastRefIdxOfWindow+1)
				}
				if lastRefIdxOfWindow != lastRefIdxOfWindow_V2 {
					fmt.Printf("lastOfWindow V1 and V2 don't agree\n")
				}
				//fmt.Printf("BeforeWindow V1 and V2 agree? %v, %v, %v\n", refIdxBeforeWindow == refIdxBeforeWindow_V2, refIdxBeforeWindow, refIdxBeforeWindow_V2)
				//fmt.Printf("lastOfWindow V1 and V2 agree? %v, %v, %v\n", lastRefIdxOfWindow == lastRefIdxOfWindow_V2, lastRefIdxOfWindow, lastRefIdxOfWindow_V2)
				// Conclusion about V2. Only last test has don't agree. Is correct window 2 12? NO I think V2 is correct

				// an option/flag can tell us not to print if there are Ns in the firstQuery or secondQuery
				if !s.RemoveN || totalNs == 0 {
					if s.LongOutput {
						percentDiverged = 100 * (float64(totalSubst+totalGaps) / float64(s.WindowSize))
						if totalSubst+totalGaps > s.WindowSize {
							log.Fatalf("Error: total number of mutations exceeds windowSize. This may or may not be a bug, but your sequence has deviated from our use case.\n")
						}
						rawPValue = scorePValueCache[totalSubst+totalGaps]
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%s\t%e\t%e\n", s.RefChromName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, s.RefChromName, refIdxBeforeWindow+1, totalSubst+totalGaps, "+", percentDiverged, rawPValue)
						exception.PanicOnErr(err)
					} else {
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\n", s.RefChromName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, s.RefChromName, refIdxBeforeWindow+1, totalSubst+totalGaps)
						// TODO: change back to code above after debugging
						//_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%d\n", s.RefChromName, refIdxBeforeWindow+1, lastRefIdxOfWindow+1, s.RefChromName, refIdxBeforeWindow+1, totalSubst+totalGaps, alnIdxBeforeWindow)
						exception.PanicOnErr(err)
					}
				}
			} else {
				//fmt.Printf("Won't report this window. alnIdxBeforeWindow+1 or lastAlnIdxOfWindow is gap in only reference\n")
			} // TODO: remove else after debugging
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
