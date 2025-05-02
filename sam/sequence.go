package sam

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// SamBedToBases takes a sam record and bed record and returns the sam sequence that corresponds to the reference coordinates
// specified by the bed record. Log-fatals if the bed record is not completely within the sam record.
func SamBedToBases(s Sam, b bed.Bed) []dna.Base {
	var idxStart int
	if !within(b, s) {
		log.Fatalf("Error: intervals don't overlap")
	}
	samPos := s.GetChromStart()
	for j := range s.Cigar {
		if s.Cigar[j].Op == 'H' {
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
			break
		}
	}
	idxEnd := idxStart + (b.ChromEnd - b.ChromStart)
	return s.Seq[idxStart:idxEnd]
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
