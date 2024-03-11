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
func incrementWindowEdge(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, alnIdxOrig int) (alnIdx, gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch int) {
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
	// 3. gapOpenedSecondQuery: did we open a gap in the secondQuery sequence when moving the window edge?
	if ((alnIdxOrig != -1 && !pDna.IsGap(secondQuery[alnIdxOrig])) || alnIdxOrig == -1) && pDna.IsGap(secondQuery[alnIdx]) {
		gapOpenedSecondQuery++
	}
	// 4. gapClosedSecondQuery: id we close a gap in the secondQuery when moving the window edge?
	if alnIdxOrig != -1 && pDna.IsGap(secondQuery[alnIdxOrig]) && !pDna.IsGap(secondQuery[alnIdx]) {
		//if seqTwo[alnIdx] == dna.Gap && (alnIdx+1==len(seqOne) || seqTwo[alnIdx+1] != dna.Gap) {
		gapClosedSecondQuery++
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

// speedyWindowDifference is a helper function of faFindFast that calculates the divergence between two input sequences (firstQuery, secondQuery sequences) using a sliding window, and then reports the divergence in terms of reference positions (positions in the reference sequence).
// optional arguments longOutput and divergenceRate allow the user to report a -log10pValue corresponding to the p value of observing a level of divergence for a given
// window under a null binomial model of neutral evolution.
func speedyWindowDifference(reference []pDna.Float32Base, firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, s Settings) {
	var alnIdxBeforeWindow, lastAlnIdxOfWindow int = -1, -1                                                                                 // these are the two edges of the sliding window in "alignment positions"
	var alnIdxBeforeWindowForRef int = -1                                                                                                   // this is to generate the alnIdx intermediate for translating firstQueryIdxBeforeWindow+1 to refIdxWindowStart
	var firstQueryIdxBeforeWindow, lastFirstQueryIdxOfWindow int = -1, -1                                                                   // these are the two edges of the sliding window in "firstQuery (no gaps) positions"
	var refIdxWindowStart, lastRefIdxOfWindowPlusOne int                                                                                    // these are the two edges of the sliding window in "reference (no gaps) positions" PLUS ONE. refIdxWindowStart = refIdxBeforeWindow+1. lastRefIdxOfWindowPlusOne = lastRefIdxOfWindow+1.
	var totalGaps, totalNs int                                                                                                              // this is the data we need to keep track of that describes the current window
	var gapOpenCloseFirstQuery, gapOpenedSecondQuery, gapClosedSecondQuery, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch int // ints we will get back when moving the window one ref base.
	var err error
	var percentDiverged, rawPValue float64
	var prevRefIdxWindowStart, prevAlnIdxBeforeWindowForRefPlusOne, prevLastRefIdxOfWindowPlusOne, prevLastAlnIdxOfWindowPlusOne int = 0, 0, 0, 0
	var windowDist []float64
	var windowDistMean float64
	var windowDistDiv int
	var windowDotSubst int // replaces totalSubst
	var windowDot []float64
	var windowDotMean float64

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
		lastAlnIdxOfWindow, gapOpenCloseFirstQuery, gapOpenedSecondQuery, _, numFirstQueryNs, numSecondQueryNsGap, numSecondQueryNsMatch = incrementWindowEdge(firstQuery, secondQuery, lastAlnIdxOfWindow)
		lastFirstQueryIdxOfWindow++
		totalGaps += gapOpenCloseFirstQuery + gapOpenedSecondQuery
		totalNs += numFirstQueryNs + numSecondQueryNsGap + numSecondQueryNsMatch
		// TODO: removed numSubst variable. Calculate totalSubst a different way

		// usually increment the baseBeforeWindow,
		// but not at the beginning when we have not yet incremented the end enough to have a full "windowSize" of bases in the window
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow > s.WindowSize {
			alnIdxBeforeWindow, _, _, _, numFirstQueryNs, _, numSecondQueryNsMatch = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow)
			alnIdxBeforeWindowForRef = updateAlnIdxBeforeWindow(firstQuery, alnIdxBeforeWindow)
			firstQueryIdxBeforeWindow++
			totalNs -= numFirstQueryNs + numSecondQueryNsMatch
		}

		// the trailing window needs to "look ahead" to see what happens before the next firstQuery base
		if lastFirstQueryIdxOfWindow-firstQueryIdxBeforeWindow == s.WindowSize {
			_, gapOpenCloseFirstQuery, _, gapClosedSecondQuery, _, numSecondQueryNsGap, _ = incrementWindowEdge(firstQuery, secondQuery, alnIdxBeforeWindow)
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

				// TODO: now that window edges are found, calculate distance score for each position and the average of the window
				// TODO: instead, calculate divergence score, both from dist and from dot
				// window edges in AlnPos [closed, open) for bed: alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow+1
				// window edges in RefPos [closed, open) for bed: refIdxWindowStart, lastRefIdxOfWindowPlusOne
				// window edges in AlnPos [closed, closed] for distance scores: alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow
				windowDist, windowDistMean = distWindow(firstQuery, secondQuery, s.WindowSize, alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow)
				baseDistToDivThreshold := 0.8
				baseDotToSubstThreshold := 0.8
				windowDistDiv = windowDistToDiv(firstQuery, secondQuery, alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow, baseDistToDivThreshold)
				windowDotSubst = windowDotToSubst(firstQuery, secondQuery, alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow, baseDotToSubstThreshold)
				windowDot, windowDotMean = dotWindow(firstQuery, secondQuery, s.WindowSize, alnIdxBeforeWindowForRef+1, lastAlnIdxOfWindow)
				// TODO: make baseDistToDivThreshold and baseDotToSubstThreshold user-input variables, from pfaFindFst, feed into speedyWindowDifference function in efficient.go

				// print output
				// an option/flag can tell us not to print if there are Ns in the firstQuery or secondQuery
				if !s.RemoveN || totalNs == 0 {
					if s.LongOutput && !s.OutputAlnPos {
						percentDiverged = 100 * (float64(windowDotSubst+totalGaps) / float64(s.WindowSize))
						//percentDiverged = 100 * (float64(totalSubst+totalGaps) / float64(s.WindowSize))
						//if totalSubst+totalGaps > s.WindowSize {
						if windowDotSubst+totalGaps > s.WindowSize {
							log.Fatalf("Error: total number of mutations exceeds windowSize. This may or may not be a bug, but your sequence has deviated from our use case.\n")
						}
						rawPValue = scorePValueCache[windowDotSubst+totalGaps]
						//rawPValue = scorePValueCache[totalSubst+totalGaps]
						// in the below output, windowDotSubst+totalGaps replaces totalSubst+totalGaps
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%s\t%e\t%e\t%d\t%e\t%e\t%d\t%e\t%e\n", s.RefChromName, refIdxWindowStart, lastRefIdxOfWindowPlusOne, s.RefChromName, refIdxWindowStart, windowDotSubst+totalGaps, "+", percentDiverged, rawPValue, windowDotSubst, windowDotMean, windowDot, windowDistDiv, windowDistMean, windowDist)
						exception.PanicOnErr(err)
					} else if !s.LongOutput && s.OutputAlnPos {
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%d\t%e\t%e\t%d\t%e\t%e\t%d\t%e\t%e\n", s.RefChromName, refIdxWindowStart, lastRefIdxOfWindowPlusOne, s.RefChromName, refIdxWindowStart, windowDotSubst+totalGaps, alnIdxBeforeWindow+1, windowDotSubst, windowDotMean, windowDot, windowDistDiv, windowDistMean, windowDist)
						exception.PanicOnErr(err)
					} else if s.LongOutput && s.OutputAlnPos {
						percentDiverged = 100 * (float64(windowDotSubst+totalGaps) / float64(s.WindowSize))
						//percentDiverged = 100 * (float64(totalSubst+totalGaps) / float64(s.WindowSize))
						//if totalSubst+totalGaps > s.WindowSize {
						// TODO: uncomment "total number of mutations exceeds windowSize" error reporting after debugging
						/*
								if windowDotSubst+totalGaps > s.WindowSize {
								log.Fatalf("Error: total number of mutations exceeds windowSize. This may or may not be a bug, but your sequence has deviated from our use case.\n")
							}
						*/
						rawPValue = scorePValueCache[windowDotSubst+totalGaps]
						//rawPValue = scorePValueCache[totalSubst+totalGaps]
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%s\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%e\t%e\n", s.RefChromName, refIdxWindowStart, lastRefIdxOfWindowPlusOne, s.RefChromName, refIdxWindowStart, windowDotSubst+totalGaps, "+", percentDiverged, rawPValue, alnIdxBeforeWindow+1, windowDotSubst, windowDotMean, windowDot, windowDistDiv, windowDistMean, windowDist)
						exception.PanicOnErr(err)
					} else {
						_, err = fmt.Fprintf(file, "%s\t%d\t%d\t%s_%d\t%d\t%d\t%e\t%e\t%d\t%e\t%e\n", s.RefChromName, refIdxWindowStart, lastRefIdxOfWindowPlusOne, s.RefChromName, refIdxWindowStart, windowDotSubst+totalGaps, windowDotSubst, windowDotMean, windowDot, windowDistDiv, windowDistMean, windowDist)
						exception.PanicOnErr(err)
					}
				}
			}
		}
	}

	// close outFile
	err = file.Close()
	exception.PanicOnErr(err)
}

func distWindow(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, windowSize int, windowStart int, windowEnd int) (windowDist []float64, windowDistMean float64) {
	var baseDist float64
	var windowDistTotal = 0.0

	for i := windowStart; i <= windowEnd; i++ {
		baseDist = pDna.Dist(firstQuery[i], secondQuery[i])
		if math.IsNaN(baseDist) {
			log.Fatalf("Error: distWindow NaN. The 2 bases are: %v, %v\n", firstQuery[i], secondQuery[i])
		}
		windowDist = append(windowDist, baseDist) // appending each baseDist instead of defining windowDist size first
		windowDistTotal += baseDist
	}

	windowDistMean = windowDistTotal / float64(windowSize)

	return windowDist, windowDistMean
}

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

func dotWindow(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, windowSize int, windowStart int, windowEnd int) (windowDot []float64, windowDotMean float64) {
	var baseDot float64
	var windowDotTotal = 0.0

	for i := windowStart; i <= windowEnd; i++ {
		baseDot = pDna.DotSubstProb(firstQuery[i], secondQuery[i])
		if math.IsNaN(baseDot) {
			log.Fatalf("Error: dotWindow NaN. The 2 bases are: %v, %v. The position is: %v\n", firstQuery[i], secondQuery[i], i)
		}
		windowDot = append(windowDot, baseDot) // appending each baseDist instead of defining windowDist size first
		windowDotTotal += baseDot
	}

	windowDotMean = windowDotTotal / float64(windowSize)

	return windowDot, windowDotMean
}

func windowDotToSubst(firstQuery []pDna.Float32Base, secondQuery []pDna.Float32Base, windowStart int, windowEnd int, baseDotToSubstThreshold float64) int {
	var baseDot float64
	var windowDotSubst = 0

	for i := windowStart; i <= windowEnd; i++ {
		if pDna.IsGap(firstQuery[i]) || pDna.IsGap(secondQuery[i]) { // do not count gaps
			continue
		} else {
			baseDot = pDna.DotSubstProb(firstQuery[i], secondQuery[i])
			if baseDot > baseDotToSubstThreshold {
				windowDotSubst++ // baseDistDiv is true, add 1 to the windowDistDiv count
			}
		}
	}

	return windowDotSubst
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
