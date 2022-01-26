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
func IlluminaSam(refName string, ref []dna.Base, numPairs int) ([]sam.Sam, sam.Header) {
	return simulatePairedSam(refName, ref, numPairs, 150, 50, 50)
}

func simulatePairedSam(refName string, ref []dna.Base, numPairs, readLen, avgInsertSize int, avgInsertSizeStdDev float64) ([]sam.Sam, sam.Header) {
	reads := make([]sam.Sam, numPairs*2)
	var header sam.Header

	var insertSize, midpoint, startFor, startRev, endFor, endRev int
	for i := 0; i < len(reads); i += 2 {
		insertSize = int(numbers.SampleInverseNormal(float64(avgInsertSize), avgInsertSizeStdDev))
		midpoint = numbers.RandIntInRange(0, len(ref))
		startFor = midpoint - (readLen + (insertSize / 2))
		endFor = startFor + readLen
		startRev = midpoint + (insertSize / 2)
		endRev = startRev + readLen

		reads[i] = generateSamReadNoFlag(fmt.Sprintf("Read:%dF", i), refName, ref, startFor, endFor)
		reads[i+1] = generateSamReadNoFlag(fmt.Sprintf("Read:%dR", i), refName, ref, startRev, endRev)
		addPairedFlags(&reads[i], &reads[i+1])
		reads[i].RNext = reads[i+1].RName
		reads[i].PNext = reads[i+1].Pos
		reads[i+1].RNext = reads[i].RName
		reads[i+1].PNext = reads[i].Pos
	}

	return reads, header
}

// soft clips sequence that is off template.
// does not generate Flag, RNext, or PNext
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

	s.MapQ = 40
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
		s.Cigar = append(s.Cigar, &cigar.Cigar{RunLength: realSeqStartIdx - start, Op: 'S', Sequence: nil})
	}
	s.Cigar = append(s.Cigar, &cigar.Cigar{RunLength: realSeqEndIdx - realSeqStartIdx, Op: 'M', Sequence: nil})
	if realSeqEndIdx < end {
		s.Cigar = append(s.Cigar, &cigar.Cigar{RunLength: end - realSeqEndIdx, Op: 'S', Sequence: nil})
	}

	return s
}

func addPairedFlags(f, r *sam.Sam) {
	var fIsRevComp bool = rand.Float64() > 0.5
	f.Flag += 1 + 40
	r.Flag += 1 + 80
	switch {
	case f.Cigar != nil && r.Cigar != nil: // both mapped
		f.Flag += 2
		r.Flag += 2
		if fIsRevComp {
			f.Flag += 10
			r.Flag += 20
		} else {
			f.Flag += 20
			r.Flag += 10
		}

	case f.Cigar == nil && r.Cigar == nil: // both unmapped
		f.Flag += 4 + 8
		r.Flag += 4 + 8

	case f.Cigar != nil && r.Cigar == nil: // f mapped r unmapped
		f.Flag += 8
		if fIsRevComp {
			f.Flag += 10
			r.Flag += 20
		}

	case f.Cigar == nil && r.Cigar != nil: // f unmapped r mapped
		r.Flag += 8
		if !fIsRevComp {
			f.Flag += 20
			r.Flag += 10
		}
	}
}
