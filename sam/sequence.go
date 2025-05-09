package sam

import (
	"errors"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// SamBedToBases takes a sam record and bed record and returns the sam sequence that corresponds to the reference coordinates
// specified by the bed record. Log-fatals if the bed record is not completely within the sam record. Currently, SamBedToBases cannot deal with any
// insertions, deletions, or introns within the specified bed region. It will return an empty slice of dna.Base and an error if it encounters this issue.
func SamBedToBases(s Sam, b bed.Bed) ([]dna.Base, error) {
	var idxStart, idxEnd int
	if !within(b, s) {
		log.Fatalf("Error in SamBedToBases: The bed interval isn't completely within the sam inerval")
	}
	samPos := s.GetChromStart()

	for j := range s.Cigar {
		if s.Cigar[j].Op == cigar.HardClip {
			continue
		}
		if cigar.ConsumesQuery(s.Cigar[j].Op) {
			idxStart += s.Cigar[j].RunLength
		}
		if cigar.ConsumesReference(s.Cigar[j].Op) {
			samPos += s.Cigar[j].RunLength
		}
		if samPos >= b.ChromStart {
			idxStart = idxStart - (samPos - b.ChromStart)
			idxEnd = idxStart + (b.ChromEnd - b.ChromStart)
			if s.Cigar[j].Op != cigar.Match || samPos < idxEnd {
				return []dna.Base{}, errors.New("error: bed interval is not in a region of alignment match and cannot be easily converted")
			}
			break
		}
	}

	return s.Seq[idxStart:idxEnd], nil
}

// within returns true if alpha falls completely within or is equal to beta, otherwise, returns false
func within(alpha bed.Bed, beta Sam) bool {
	if alpha.GetChrom() != beta.GetChrom() {
		return false
	}
	if alpha.GetChromStart() >= beta.GetChromStart() && alpha.GetChromEnd() <= beta.GetChromEnd() {
		return true
	}
	return false
}
