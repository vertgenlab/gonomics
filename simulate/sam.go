package simulate

import (
	"fmt"
	"math/rand"
	"strings"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
)

// IlluminaPairedSam generates a pair of sam reads randomly distributed across the input ref.
func IlluminaPairedSam(refName string, ref []dna.Base, numPairs, readLen, avgFragmentSize int, avgFragmentStdDev float64, flatErrorRate float64, binomialAlias numbers.BinomialAlias, out *fileio.EasyWriter, bw *sam.BamWriter, bamOutput bool) {
	var fragmentSize, midpoint, startFor, startRev, endFor, endRev int
	var currFor, currRev sam.Sam
	for i := 0; i < numPairs; i++ {
		fragmentSize = numbers.Max(readLen, int(numbers.SampleInverseNormal(float64(avgFragmentSize), avgFragmentStdDev)))
		midpoint = numbers.RandIntInRange(0, len(ref))
		startFor = midpoint - (fragmentSize / 2)
		endFor = startFor + readLen
		endRev = midpoint + (fragmentSize / 2)
		startRev = endRev - readLen

		currFor = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i), refName, ref, startFor, endFor, flatErrorRate, binomialAlias)
		currRev = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i), refName, ref, startRev, endRev, flatErrorRate, binomialAlias)
		if currFor.Cigar == nil && currRev.Cigar == nil {
			i -= 1 // retry
			continue
		}
		addPairedFlags(&currFor, &currRev)
		if currFor.Cigar != nil && currRev.Cigar != nil {
			currFor.RNext = "="
			currRev.RNext = "="
		} else {
			currFor.RNext = currRev.RName
			currRev.RNext = currFor.RName
		}

		currFor.PNext = currRev.Pos
		currRev.PNext = currFor.Pos
		if bamOutput {
			sam.WriteToBamFileHandle(bw, currFor, 0)
			sam.WriteToBamFileHandle(bw, currRev, 0)
		} else {
			sam.WriteToFileHandle(out, currFor)
			sam.WriteToFileHandle(out, currRev)
		}
	}
}

// generateSamReadNoFlag generates a sam record for the input position.
// Soft clips sequence that is off template and does not generate Flag, RNext, or PNext.
func generateSamReadNoFlag(readName string, refName string, ref []dna.Base, start, end int, flatErrorRate float64, alias numbers.BinomialAlias) sam.Sam {
	var currSam sam.Sam
	currSam.QName = readName
	currSam.Seq = make([]dna.Base, end-start)
	// generate qual
	var bldr strings.Builder
	for range currSam.Seq {
		bldr.WriteRune(rune(numbers.RandIntInRange(30, 40) + 33)) // high quality seq + ascii offset
	}
	currSam.Qual = bldr.String()

	// check if unmapped
	if end < 0 || start > len(ref) {
		currSam.RName = "*"
		for i := range currSam.Seq {
			currSam.Seq[i] = dna.Base(numbers.RandIntInRange(0, 4))
		}
		return currSam
	}

	currSam.MapQ = uint8(numbers.RandIntInRange(30, 40))
	currSam.RName = refName

	// generate random seq if off template
	var realSeqStartIdx, realSeqEndIdx int
	for realSeqStartIdx = start; realSeqStartIdx < 0; realSeqStartIdx++ {
		currSam.Seq[realSeqStartIdx-start] = dna.Base(numbers.RandIntInRange(0, 4))
	}
	for realSeqEndIdx = end; realSeqEndIdx > len(ref); realSeqEndIdx-- {
		currSam.Seq[len(currSam.Seq)-(1+(end-realSeqEndIdx))] = dna.Base(numbers.RandIntInRange(0, 4))
	}
	copy(currSam.Seq[realSeqStartIdx-start:len(currSam.Seq)-(end-realSeqEndIdx)], ref[realSeqStartIdx:realSeqEndIdx])

	// now we will simulate sequencing/PCR error with a flat error rate
	if flatErrorRate > 0 {
		currSam = sequencingError(currSam, start, end, alias)
	}

	// generate other values
	currSam.Pos = uint32(realSeqStartIdx) + 1
	currSam.TLen = int32(realSeqEndIdx - realSeqStartIdx)

	// assemble cigar
	if realSeqStartIdx > start {
		currSam.Cigar = append(currSam.Cigar, cigar.Cigar{RunLength: realSeqStartIdx - start, Op: 'S'})
	}
	currSam.Cigar = append(currSam.Cigar, cigar.Cigar{RunLength: realSeqEndIdx - realSeqStartIdx, Op: 'M'})
	if realSeqEndIdx < end {
		currSam.Cigar = append(currSam.Cigar, cigar.Cigar{RunLength: end - realSeqEndIdx, Op: 'S'})
	}

	return currSam
}

// addPairedFlags adds the flag for a pair of sam records.
func addPairedFlags(f, r *sam.Sam) {
	var fIsRevComp bool = rand.Float64() > 0.5
	if fIsRevComp {
		*f, *r = *r, *f // so that the reads always point towards one another
	}
	f.Flag += 1 + 64
	r.Flag += 1 + 128
	switch {
	case f.Cigar != nil && r.Cigar != nil: // both mapped
		f.Flag += 2
		r.Flag += 2
		if fIsRevComp {
			f.Flag += 16
			r.Flag += 32
		} else {
			f.Flag += 32
			r.Flag += 16
		}

	case f.Cigar == nil && r.Cigar == nil: // both unmapped
		f.Flag += 4 + 8
		r.Flag += 4 + 8

	case f.Cigar != nil && r.Cigar == nil: // f mapped r unmapped
		f.Flag += 8
		r.Flag += 4
		if fIsRevComp {
			f.Flag += 16
			r.Flag += 32
		}

	case f.Cigar == nil && r.Cigar != nil: // f unmapped r mapped
		f.Flag += 4
		r.Flag += 8
		if !fIsRevComp {
			f.Flag += 32
			r.Flag += 16
		}
	}
}

// sequencingError takes in a sam.Sam record, start and end positions, and a BinomialAlias to
// edit bases randomly across the read to simulate PCR/sequencing error. Errors are made with no
// positional dependence and with a flat error spectrum.
func sequencingError(currSam sam.Sam, start int, end int, alias numbers.BinomialAlias) sam.Sam {
	numFlatErrors := numbers.RandBinomial(alias)         // sample a binomial distribution to get the number of sequencing errors
	mutatedPositions := make(map[int]int, numFlatErrors) // store positions we've mutated so that we can sample without replacement
	var foundInMap bool
	var currRandInt int
	var currError int = 0

	for currError < numFlatErrors {
		currRandInt = numbers.RandIntInRange(0, end-start) // sample a base on the read
		if _, foundInMap = mutatedPositions[currRandInt]; !foundInMap {
			mutatedPositions[currRandInt] = 1
			currSam.Seq[currRandInt] = changeBase(currSam.Seq[currRandInt])
			currError++
		}
	}
	return currSam
}
