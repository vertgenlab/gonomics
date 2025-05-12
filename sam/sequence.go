package sam

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// SamBedToBases takes a sam record and bed record and returns the sam sequence that corresponds to the reference coordinates
// specified by the bed record. Log-fatals if the bed record is not completely within the sam record. Currently, SamBedToBases only returns bases that match the reference,
// (that have an M annotation in the cigar). If other cigar annotations are present, the output will be truncated to the matching portion. Additionally, an "annotation" string corresponding to the cigar
// values behind the sequence in the output, one letter in the annotation string per base. At a minimum the string will have the same length as the bed region (using 'D' or 'N') if there are introns or deletion that truncate the sequence output.
// The annotation string can also be longer than the bed region in the case on a sequence insertion in the Sam record. In that case 'I' will be used in the annotation string.
func SamBedToBases(s Sam, b bed.Bed) (ans []dna.Base, annotation string) {
	var idx, idxStart, idxEnd int
	var withinBed, secondIter bool = false, false
	if !within(b, s) {
		log.Fatalf("Error in SamBedToBases: The bed interval isn't completely within the sam inerval")
	}
	samPos := s.GetChromStart()

	for j := range s.Cigar {
		if s.Cigar[j].Op == cigar.HardClip {
			continue
		}
		//update sam sequence index and sam sequence position
		if cigar.ConsumesQuery(s.Cigar[j].Op) {
			idx += s.Cigar[j].RunLength
		}
		if cigar.ConsumesReference(s.Cigar[j].Op) {
			samPos += s.Cigar[j].RunLength
		}

		//for building the annotation string, we need to know if the loop is on an iteration after we have been within the bed region
		if samPos >= b.ChromStart && withinBed {
			secondIter = true
		}

		//if sam position is above the bed start value, we need to set the start index of the sequence we want to grab
		if samPos >= b.ChromStart && !withinBed {
			withinBed = true
			if s.Cigar[j].Op == cigar.Match {
				idxStart = idx - (samPos - b.ChromStart)
			} else {
				idxStart = idx
			}

			//update the annotation string with different rules based on if this will be the last loop in the interation
			if samPos >= b.ChromEnd {
				annotation = buildAnnoString(annotation, string(s.Cigar[j].Op), b.ChromEnd-b.ChromStart)
			} else {
				annotation = buildAnnoString(annotation, string(s.Cigar[j].Op), samPos-b.ChromStart)
			}

		}
		//if sam position in above the bed end value, we need to set the end index of the sequence we want to grab
		if samPos >= b.ChromEnd {
			if s.Cigar[j].Op == cigar.Match {
				idxEnd = idx - (samPos - b.ChromEnd)
			} else {
				idxEnd = idx
			}

			//update the annotation string if this isn't the only iteration that overlaps the bed region (if it is, the annotation string has already been build above)
			if secondIter {
				annotation = buildAnnoString(annotation, string(s.Cigar[j].Op), s.Cigar[j].RunLength-(samPos-b.ChromEnd))
			}
			break
		}

		//if we are in the bed region for the 2+ time, but we still haven't reached the end of the bed region, we need to update the annotation string for the current cigar operation and run length
		if secondIter {
			annotation = buildAnnoString(annotation, string(s.Cigar[j].Op), s.Cigar[j].RunLength)
		}
	}

	return s.Seq[idxStart:idxEnd], annotation
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

func buildAnnoString(anno string, op string, count int) string {
	for i := 0; i < count; i++ {
		anno += op
	}
	return anno
}
