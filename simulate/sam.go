package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"math/rand"
	"strings"
)

// IlluminaSam generates pseudorandom paired reads evenly distributed along the input ref sequence using typical illumina parameters.
func IlluminaSam(refName string, ref []dna.Base, numPairs int) []sam.Sam {
	return simulatePairedSam(refName, ref, numPairs, 150, 50, 50)
}

// simulatePairedSam generates a pair of sam reads randomly distributed across the input ref.
func simulatePairedSam(refName string, ref []dna.Base, numPairs, readLen, avgInsertSize int, avgInsertSizeStdDev float64) []sam.Sam {
	reads := make([]sam.Sam, numPairs*2)

	var insertSize, midpoint, startFor, startRev, endFor, endRev int
	for i := 0; i < len(reads); i += 2 {
		insertSize = int(numbers.SampleInverseNormal(float64(avgInsertSize), avgInsertSizeStdDev))
		midpoint = numbers.RandIntInRange(0, len(ref)) // tapered coverage at ends of contig
		//midpoint = numbers.RandIntInRange(-2*readLen, len(ref)+(2*readLen)) // even coverage across entire contig
		startFor = midpoint - (readLen + (insertSize / 2))
		endFor = startFor + readLen
		startRev = midpoint + (insertSize / 2)
		endRev = startRev + readLen

		reads[i] = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i/2), refName, ref, startFor, endFor)
		reads[i+1] = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i/2), refName, ref, startRev, endRev)
		if reads[i].Cigar == nil && reads[i+1].Cigar == nil {
			i -= 2 // retry
			continue
		}
		addPairedFlags(&reads[i], &reads[i+1])
		if reads[i].Cigar != nil && reads[i+1].Cigar != nil {
			reads[i].RNext = "="
			reads[i+1].RNext = "="
		} else {
			reads[i].RNext = reads[i+1].RName
			reads[i+1].RNext = reads[i].RName
		}

		reads[i].PNext = reads[i+1].Pos
		reads[i+1].PNext = reads[i].Pos
	}
	return reads
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
