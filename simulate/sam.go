package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"math/rand"
	"strings"
)

// simulatePairedSam generates a pair of sam reads randomly distributed across the input ref.
func IlluminaPairedSam(refName string, ref []dna.Base, numPairs, readLen, avgFragmentSize int, avgFragmentStdDev float64, out *fileio.EasyWriter, bw *sam.BamWriter, bamOutput bool) {
	var fragmentSize, midpoint, startFor, startRev, endFor, endRev int
	var currFor, currRev sam.Sam
	for i := 0; i < numPairs; i++ {
		fragmentSize = int(numbers.SampleInverseNormal(float64(avgFragmentSize), avgFragmentStdDev))
		midpoint = numbers.RandIntInRange(0, len(ref))
		startFor = midpoint - (fragmentSize / 2)
		endFor = startFor + readLen
		endRev = midpoint + (fragmentSize / 2)
		startRev = endRev - readLen

		currFor = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i), refName, ref, startFor, endFor)
		currRev = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i), refName, ref, startRev, endRev)
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
func generateSamReadNoFlag(readName string, refName string, ref []dna.Base, start, end int) sam.Sam {
	var s sam.Sam
	s.QName = readName
	s.Seq = make([]dna.Base, end-start)
	// generate qual
	var bldr strings.Builder
	for range s.Seq {
		bldr.WriteRune(rune(numbers.RandIntInRange(30, 40) + 33)) // high quality seq + ascii offset
	}
	s.Qual = bldr.String()

	// check if unmapped
	if end < 0 || start > len(ref) {
		s.RName = "*"
		for i := range s.Seq {
			s.Seq[i] = dna.Base(numbers.RandIntInRange(0, 4))
		}
		return s
	}

	s.MapQ = uint8(numbers.RandIntInRange(30, 40))
	s.RName = refName

	// generate random seq if off template
	var realSeqStartIdx, realSeqEndIdx int
	for realSeqStartIdx = start; realSeqStartIdx < 0; realSeqStartIdx++ {
		s.Seq[realSeqStartIdx-start] = dna.Base(numbers.RandIntInRange(0, 4))
	}
	for realSeqEndIdx = end; realSeqEndIdx > len(ref); realSeqEndIdx-- {
		s.Seq[len(s.Seq)-(1+(end-realSeqEndIdx))] = dna.Base(numbers.RandIntInRange(0, 4))
	}
	copy(s.Seq[realSeqStartIdx-start:len(s.Seq)-(end-realSeqEndIdx)], ref[realSeqStartIdx:realSeqEndIdx])

	// generate other values
	s.Pos = uint32(realSeqStartIdx) + 1
	s.TLen = int32(realSeqEndIdx - realSeqStartIdx)

	// assemble cigar
	if realSeqStartIdx > start {
		s.Cigar = append(s.Cigar, cigar.Cigar{RunLength: realSeqStartIdx - start, Op: 'S'})
	}
	s.Cigar = append(s.Cigar, cigar.Cigar{RunLength: realSeqEndIdx - realSeqStartIdx, Op: 'M'})
	if realSeqEndIdx < end {
		s.Cigar = append(s.Cigar, cigar.Cigar{RunLength: end - realSeqEndIdx, Op: 'S'})
	}

	return s
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
