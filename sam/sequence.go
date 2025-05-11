package sam

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// SamBedToBases takes a sam record and bed record and returns the sam sequence that corresponds to the reference coordinates
// specified by the bed record. Log-fatals if the bed record is not completely within the sam record. Currently, SamBedToBases only returns bases that match the reference,
// (that have an M annotation in the cigar). If other cigar annotations are present, the output will be truncated to the matching portion
func SamBedToBases(s Sam, b bed.Bed) ([]dna.Base, error) {
	var idx, idxStart, idxEnd int
	var withinBed bool = false
	if !within(b, s) {
		log.Fatalf("Error in SamBedToBases: The bed interval isn't completely within the sam inerval")
	}
	samPos := s.GetChromStart()

	for j := range s.Cigar {
		if s.Cigar[j].Op == cigar.HardClip {
			continue
		}
		if cigar.ConsumesQuery(s.Cigar[j].Op) {
			idx += s.Cigar[j].RunLength
		}
		if cigar.ConsumesReference(s.Cigar[j].Op) {
			samPos += s.Cigar[j].RunLength
		}
		if samPos >= b.ChromStart && !withinBed {
			withinBed = true
			if s.Cigar[j].Op == cigar.Match {
				idxStart = idx - (samPos - b.ChromStart)
			} else {
				idxStart = idx
			}
		}
		if samPos >= b.ChromEnd {
			if s.Cigar[j].Op == cigar.Match {
				idxEnd = idx - (samPos - b.ChromEnd)
			} else {
				idxEnd = idx
			}
			break
		}
	}
	//fmt.Println(idx)
	//fmt.Println(samPos)
	//fmt.Println(b.ChromStart, b.ChromEnd)
	//fmt.Println(idxStart, idxEnd)
	if !withinBed {
		idxStart = idx
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
