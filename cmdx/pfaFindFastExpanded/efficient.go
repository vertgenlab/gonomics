package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta/pFasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/logspace"
	"log"
	"math"
	"strings"
)

// incrementWindowEdge
// inputs: 2 aligned pFasta sequences (firstQuery, secondQuery), a current alignment position index.
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
//  9. numConfident: the number of secondQuery positions where the most likely base >= threshold, within numSubst (<= numSubst), float32 threshold
func incrementWindowEdge(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, alnIdxOrig int, baseDotToSubstThreshold float64, confidentThreshold float32) (alnIdx, gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst, numConfident int) {
	alnIdx = alnIdxOrig
	// at initialization, lastAlnIdxOfWindow == -1, which is assigned to alnIdxOrig when incrementWindowEdge is called, so alnIdx = alnIdxOrig == -1. For loop will start at alnIdx++ which is 0, the 1st position in the alignment

	// increment alnIdx to the next non-gap position in the firstQuery sequence
	// for loop start condition: start at alnIdx++, which is the position after the current alignment position
	// for loop end condition: ends when we reach the end of the firstQuery sequence, or when we reach a non-gap position in the firstQuery sequence, whichever happens first
	// the for loop's start and end conditions also mean that the for loop is never executed, unless there is firstQuery gap before we reach the end of firstQuery
	// for loop incrementation: alnIdx++
	for alnIdx++; alnIdx < len(firstQuery) && pDna.IsGap(firstQuery[alnIdx]); alnIdx++ {
		// keep track of changes between the original and next non-gap positions in the firstQuery sequence
		// 6. secondQuery
		if pDna.IsN(secondQuery[alnIdx]) {
			numSecondQueryNsGap++
		}
		// 2. gapOpenClosedFirstQuery
		if !pDna.IsGap(secondQuery[alnIdx]) {
			gapOpenCloseFirstQuery = 1 // Rationale: if the for loop is executed, this means there is a gap in the firstQuery, and if secondQuery does not have a gap, then it is indeed a firstQuery gap compared to secondQuery
		}
	}

	// if we ran off the end of firstQuery when looking for the next non-gap firstQuery base, aka the for loop ended because we reached the end of the firstQuery sequence
	if alnIdx == len(firstQuery) {
		return // return all 8 named return variables
	}

	// 5. numFirstQueryNs: did we add another firstQuery N after moving the window edge by one firstQuery base
	if pDna.IsN(firstQuery[alnIdx]) {
		numFirstQueryNs++
	}
	// 7. numSecondQueryNsMatch: did we add another N to the secondQuery count of Ns after moving the window edge by one firstQuery base
	if pDna.IsN(secondQuery[alnIdx]) {
		numSecondQueryNsMatch++ // Rationale: these secondQueryNs align with non-gap firstQuery positions
	}
	// 8. numSubst: is this a substitution?
	var baseDot float64
	if !pDna.IsGap(firstQuery[alnIdx]) && !pDna.IsGap(secondQuery[alnIdx]) { // do not calculate over gaps
		baseDot = pDna.DotSubstProb(firstQuery[alnIdx], secondQuery[alnIdx]) // for non-gap position, calculate substitution probability from dot product method
		if baseDot >= baseDotToSubstThreshold {                              // if the substitution probability >= threshold
			numSubst++ // then this position is a substitution
		}
	}
	// 3. gapOpenedSecondQuery: did we open a gap in the secondQuery sequence when moving the window edge?
	if ((alnIdxOrig != -1 && !pDna.IsGap(secondQuery[alnIdxOrig])) || alnIdxOrig == -1) && pDna.IsGap(secondQuery[alnIdx]) {
		gapOpenedSecondQuery++
	}
	// 4. gapClosedSecondQuery: id we close a gap in the secondQuery when moving the window edge?
	if alnIdxOrig != -1 && pDna.IsGap(secondQuery[alnIdxOrig]) && !pDna.IsGap(secondQuery[alnIdx]) {
		//if seqTwo[alnIdx] == dna.Gap && (alnIdx+1==len(seqOne) || seqTwo[alnIdx+1] != dna.Gap) {
		gapClosedSecondQuery++
	}

	// 9. numConfident
	if pDna.IsConfident(secondQuery[alnIdx], confidentThreshold) {
		numConfident++
	}
	return
}

// updateAlnIdxBeforeWindow
// When reporting in firstQuery positions,
// alnIdxBeforeWindow corresponds to firstQueryIdxBeforeWindow, which is non-gap in firstQuery
// output bed reports chromStart = firstQueryIdxBeforeWindow+1
// However, firstQueryIdxBeforeWindow+1 may correspond not to alnIdxBeforeWindow+1
// If there are gaps in the multiFa alignment in the firstQuery,
// then firstQueryIdxBeforeWindow+1 corresponds to alnIdxBeforeWindow+(some number > 1)
// When reporting in reference positions,
// chromStart = refIdxWindowStart, which is converted from the alnIdxBeforeWindow+(some number > 1)
// in order to match firstQueryIdxBeforeWindow+1
// This function uses a while loop to find alnIdxBeforeWindow+(some number > 1), which corresponds to firstQueryIdxBeforeWindow+1
// Note that this function is only needed for chromStart, not chromEnd
// firstQuery[lastAlnIdxOfWindow+1] is allowed to be gap or exceed the firstQuery length
// This is because the bed's chromStart is closed, meaning the position firstQueryIdxBeforeWindow+1 is included in the window (for all species)
// But the bed's chromEnd is open, meaning the position lastAlnIdxOfWindow+1 is not included in the window (for all species)
// Note also that updateAlnIdxBeforeWindow is only for translating positions and does not affect divergence calculations
// BeforeWindow works better than trying to report WindowStart, because +1 is only needed once outside this function
func updateAlnIdxBeforeWindow(firstQuery []pDna.Float32Base, alnIdxOrig int) (alnIdx int) {
	alnIdx = alnIdxOrig
	for alnIdx+1 < len(firstQuery) && pDna.IsGap(firstQuery[alnIdx+1]) {
		alnIdx++
	}
	return
}

// speedyWindowDifference is a helper function of pfaFindFast that calculates the divergence between two input sequences (firstQuery, secondQuery sequences) using a sliding window, and then reports the divergence in terms of reference positions (positions in the reference sequence).
// optional arguments longOutput and divergenceRate allow the user to report a -log10pValue corresponding to the p value of observing a level of divergence for a given
// window under a null binomial model of neutral evolution.
func speedyWindowDifference(reference []pDna.Float32Base, firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, s Settings) {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1                                                                                                         // these are the two edges of the sliding window in "alignment positions"
	var alnIdxBeforeWindowForRef int = -1                                                                                                                           // this is to generate the alnIdx intermediate for translating firstQueryIdxBeforeWindow+1 to refIdxWindowStart
	var firstQueryIdxBeforeWindow, lastFirstQueryIdxOfWindow int = -1, -1                                                                                           // these are the two edges of the sliding window in "firstQuery (no gaps) positions"
	var refIdxWindowStart, lastRefIdxOfWindowPlusOne int                                                                                                            // these are the two edges of the sliding window in "reference (no gaps) positions" PLUS ONE. refIdxWindowStart = refIdxBeforeWindow+1. lastRefIdxOfWindowPlusOne = lastRefIdxOfWindow+1.
	var totalGaps, totalNs, totalSubst, totalConfident int                                                                                                          // this is the data we need to keep track of that describes the current window
	var gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst, numConfident int // ints we will get back when moving the window one ref base.
	var err error
	var percentDiverged, rawPValue float64
	var prevRefIdxWindowStart, prevAlnIdxBeforeWindowForRefPlusOne, prevLastRefIdxOfWindowPlusOne, prevLastAlnIdxOfWindowPlusOne int = 0, 0, 0, 0

	// this caches map[k] to -log10(BinomialDist(n, k, p, true)), which is the -log10 p Value.
	var scorePValueCache map[int]float64

	// create outFile
	file := fileio.EasyCreate(s.OutFile)

	// divergenceRate = math.MaxFloat64 is a reserved value that signifies that the user has not set a divergence rate. If divergenceRate != -1,
	// we initialize the scorePValueCache.
	if s.DivergenceRate != math.MaxFloat64 {
		scorePValueCache = binomialDistCacheLog10(s.WindowSize, s.DivergenceRate)
	}

	for lastAlnIdxOfWindow < len(firstQuery) { // this check could also be "!done", I am not sure what is more clear
		// we always move the lastBaseOfTheWindow (right side) and add on what we find to the counters because
		// all this stuff is now inside the current window
		lastAlnIdxOfWindow, gapOpenCloseFirstQuery, gapOpenedSecondQuery, _, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch, numSubst, numConfident = incrementWindowEdge(firstQuery, secondQuery, lastAlnIdxOfWindow, s.BaseDotToSubstThreshold, s.ConfidentThreshold)
		lastFirstQueryIdxOfWindow++
		totalGaps += gapOpenCloseFirstQuery + gapOpenedSecondQuery
		totalNs += numFirstQueryNs + numSecondQueryNsGap + numSecondQueryNsMatch
		totalSubst += numSubst
		totalConfident += numConfident

		// usually increment the baseBeforeWindow,
		// but not at the beginning when we have not yet incremented the end enough to have a full "windowSize" of bases in the window
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow > s.WindowSize {
			alnIdxBeforeWindow, _, _, _, numFirstQueryNs, _, numSecondQueryNsMatch, numSubst, numConfident = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow, s.BaseDotToSubstThreshold, s.ConfidentThreshold)
			alnIdxBeforeWindowForRef = updateAlnIdxBeforeWindow(firstQuery, alnIdxBeforeWindow)
			firstQueryIdxBeforeWindow++
			totalNs -= numFirstQueryNs + numSecondQueryNsMatch
			totalSubst -= numSubst
			totalConfident -= numConfident
		}

		// the trailing window needs to "look ahead" to see what happens before the next firstQuery base
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow == s.WindowSize {
			_, gapOpenCloseFirstQuery, _, gapClosedSecondQuery, _, numSecondQueryNsGap, _, _, _ = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow, s.BaseDotToSubstThreshold, s.ConfidentThreshold)
			totalGaps -= gapOpenCloseFirstQuery + gapClosedSecondQuery
			totalNs -= numSecondQueryNsGap
		}

		// we check to make sure we are not at the very beginning or end, where we would have partial or illegal windows
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow == s.WindowSize && lastAlnIdxOfWindow < len(firstQuery) {

			// convert position from firstQuery to reference
			// check if firstQuery chromStart (non-gap) is gap in reference. If so, don't report that window
			// check if firstQuery chromEnd-1 (last non-gap position) is gap in reference. If so, don't report that window
			if !((pDna.IsGap(reference[alnIdxBeforeWindowForRef+1]) && !pDna.IsGap(firstQuery[alnIdxBeforeWindowForRef+1])) || (pDna.IsGap(reference[lastAlnIdxOfWindow]) && !pDna.IsGap(firstQuery[lastAlnIdxOfWindow]))) {
				// call fasta.AlnPosToRefPosCounterSeq on the alnPos of chromStart and chromEnd
				// as opposed to alnIdxBeforeWindowForRef and lastAlnIdxOfWindow
				// this is because fasta.AlnPosToRefPosCounterSeq rounds position up when it encounters gap
				// in order to only scan the genome once, call fasta.AlnPosToRefPosCounterSeq on saved refStart and alnStart once refStart and alnStart both exceed 0
				if lastRefIdxOfWindowPlusOne < 0 || lastAlnIdxOfWindow+1 < 0 {
					refIdxWindowStart = pFasta.PAlnPosToRefPosCounterSeq(reference, alnIdxBeforeWindowForRef+1, 0, 0)
					lastRefIdxOfWindowPlusOne = pFasta.PAlnPosToRefPosCounterSeq(reference, lastAlnIdxOfWindow+1, refIdxWindowStart, alnIdxBeforeWindowForRef+1)
				} else {
					refIdxWindowStart = pFasta.PAlnPosToRefPosCounterSeq(reference, alnIdxBeforeWindowForRef+1, prevRefIdxWindowStart, prevAlnIdxBeforeWindowForRefPlusOne)
					lastRefIdxOfWindowPlusOne = pFasta.PAlnPosToRefPosCounterSeq(reference, lastAlnIdxOfWindow+1, prevLastRefIdxOfWindowPlusOne, prevLastAlnIdxOfWindowPlusOne)
				}

				prevRefIdxWindowStart = refIdxWindowStart
				prevAlnIdxBeforeWindowForRefPlusOne = alnIdxBeforeWindowForRef + 1
				prevLastRefIdxOfWindowPlusOne = lastRefIdxOfWindowPlusOne
				prevLastAlnIdxOfWindowPlusOne = lastAlnIdxOfWindow + 1

				// now that window edges are found, calculate divergence score (distance or dot product) for each position, the average of the window, and conversion to substitution
				// window edges in AlnPos [closed, open) for bed: alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow+1
				// window edges in RefPos [closed, open) for bed: refIdxWindowStart, lastRefIdxOfWindowPlusOne
				// window edges in AlnPos [closed, closed] for distance scores: alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow
				// print output
				// an option/flag can tell us not to print if there are Ns in the firstQuery or secondQuery
				if !s.RemoveN || totalNs == 0 {
					// use strings.Builder to construct the bed output line for each window
					var sb strings.Builder

					// common prefix: chrom, refStart, refEnd, composite name, and mutation count
					//_, err = fmt.Fprintf(&sb, "%s\t%d\t%d\t%s_%d\t%d",
					_, err = fmt.Fprintf(&sb, "%s\t%d\t%d\t%s_%d\t%d\t%d\t%d\t%d\t%d",
						s.RefChromName,
						refIdxWindowStart,
						lastRefIdxOfWindowPlusOne,
						s.RefChromName,
						refIdxWindowStart,
						totalSubst+totalGaps,
						totalSubst, //TODO: for exapnded
						totalGaps,
						lastRefIdxOfWindowPlusOne-refIdxWindowStart,
						totalConfident,
					)

					// Long output needs percentDiverged and rawPValue
					if s.LongOutput {
						percentDiverged = 100 * (float64(totalSubst+totalGaps) / float64(s.WindowSize))
						if totalSubst+totalGaps > s.WindowSize {
							log.Fatalf("Error: total number of mutations exceeds windowSize. This may or may not be a bug, but your sequence has deviated from our use case.\n")
						}
						// rawPValue only valid if cache exists; original code expects it when DivergenceRate was set
						if scorePValueCache != nil {
							rawPValue = scorePValueCache[totalSubst+totalGaps]
						} else {
							rawPValue = 0
						}
						// add strand, percent and pvalue
						// in the below output, windowDotSubst+totalGaps replaces totalSubst+totalGaps
						_, err = fmt.Fprintf(&sb, "\t%s\t%e\t%e", "+", percentDiverged, rawPValue)
					}

					// Output alignment position if requested
					if s.OutputAlnPos {
						// report aln position
						_, err = fmt.Fprintf(&sb, "\t%d", alnIdxBeforeWindow+1)
					}

					// finish line and write once
					sb.WriteByte('\n')
					_, err = fmt.Fprintf(file, sb.String())
					exception.PanicOnErr(err)
				}
			}
		}
	}

	// close outFile
	err = file.Close()
	exception.PanicOnErr(err)
}

/*
// windowDistToDiv takes 2 pDna sequences (aligned in a multi-pFa format), window start and end positions, and a user-specified float64 threshold as inputs.
// The function counts the number of divergences between the 2 pDna sequences within the window.
// The rules for counting divergences are as follows:
// Gaps are not counted as divergences.
// If the Euclidean distance between the 2 pDna base vectors at the same position in the 2 pDna sequences is greater than the threshold, then a divergence is counted.
func windowDistToDiv(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, windowStart int, windowEnd int, baseDistToDivThreshold float64) int {
	var baseDist float64
	var windowDistDiv = 0

	for i := windowStart; i <= windowEnd; i++ {
		if pDna.IsGap(firstQuery[i]) || pDna.IsGap(secondQuery[i]) { // do not count gaps
			continue
		} else {
			baseDist = pDna.Dist(firstQuery[i], secondQuery[i])
			if baseDist > baseDistToDivThreshold {
				windowDistDiv++ // baseDistDiv is true, add 1 to the windowDistDiv count
			}
		}
	}

	return windowDistDiv
}

// windowDotToSubst takes 2 pDna sequences (aligned in a multi-pFa format), window start and end positions, and a user-specified float64 threshold as inputs.
// The function counts the number of substitutions between the 2 pDna sequences within the window.
// The rules for counting substitutions are as follows:
// Gaps are not counted as substitutions.
// If the dot product distance (1 - dot product) between the 2 pDna base vectors at the same position in the 2 pDna sequences is greater than the threshold, then a substitution is counted.
func windowDotToSubst(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, windowStart int, windowEnd int, baseDotToSubstThreshold float64) int {
	var baseDot float64
	var windowDotSubst = 0

	for i := windowStart; i <= windowEnd; i++ {
		if pDna.IsGap(firstQuery[i]) || pDna.IsGap(secondQuery[i]) { // do not count gaps
			continue
		} else {
			baseDot = pDna.DotSubstProb(firstQuery[i], secondQuery[i])
			if baseDot > baseDotToSubstThreshold {
				//fmt.Printf("windowDotToSubst. firstQuery[i]: %v, secondQuery[i]: %v, baseDot: %v, i: %v\n", firstQuery[i], secondQuery[i], baseDot, i) // uncomment this line and run pfaFindFast on short sequence (not whole chr) to see pFasta calculations
				windowDotSubst++ // baseDistDiv is true, add 1 to the windowDistDiv count
			}
		}
	}

	return windowDotSubst
}
*/

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
