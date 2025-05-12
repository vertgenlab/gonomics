package sam

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

// SamBedToBases takes a sam record and bed record and returns the sam sequence that corresponds to the reference coordinates
// specified by the bed record. The program log-fatal's if the bed record is not completely within the sam record. If deletions or introns exist in the sam sequence,
// they will be annotated with gaps in the output (-). If insertions exist in the sam sequence they will be represented with lower-case bases in the output. The length
// of the output will be the same length as the bed record at a minimum, but can be longer if insertions exist in the sam sequence within the bed interval.
func SamBedToBases(s Sam, b bed.Bed) []dna.Base {
	var idx, idxStart, idxEnd int
	var withinBed bool = false
	var ans, tmp []dna.Base

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

		//if the above position/index updates occur in the loop after we have already crossed into bed region, we go in here
		if withinBed {
			//if the position/index updates take us past the end of the bed region, we need to grab the remaining sequence/gaps and exit the loop
			if samPos >= b.ChromEnd {
				idxEnd = idx - (samPos - b.ChromEnd)
				idxStart = idx - s.Cigar[j].RunLength
				if s.Cigar[j].Op == cigar.Match {
					return append(ans, s.Seq[idxStart:idxEnd]...)
				} else {
					return addGaps(ans, idxEnd-idxStart)
				}
			}

			//if we still aren't passed, we need to add the sequence/gaps pertaining to the cigar value
			if s.Cigar[j].Op == cigar.Match {
				idxStart = idx - s.Cigar[j].RunLength
				ans = append(ans, s.Seq[idxStart:idx]...)
			} else if s.Cigar[j].Op == cigar.Insertion {
				tmp = s.Seq[idx-s.Cigar[j].RunLength : idx]
				dna.AllToLower(tmp)
				ans = append(ans, tmp...)
			} else if cigar.ConsumesReference(s.Cigar[j].Op) {
				ans = addGaps(ans, s.Cigar[j].RunLength)
			}
			continue
		}

		//if sam position is above the bed start value, we need to set the start index of the sequence we want to grab
		if samPos >= b.ChromStart && !withinBed {
			withinBed = true
			//get the start sequence index of the bed region
			idxStart = idx - (samPos - b.ChromStart)

			//simplest case, grab the sequence and return
			idxStart = idx - (samPos - b.ChromStart)
			if samPos >= b.ChromEnd {
				if s.Cigar[j].Op == cigar.Match {
					idxEnd = idx - (samPos - b.ChromEnd)
					return s.Seq[idxStart:idxEnd]
				} else if s.Cigar[j].Op == cigar.Deletion || s.Cigar[j].Op == cigar.Ns {
					return addGaps(ans, b.ChromEnd-b.ChromStart)
				}
			}
			//if there are multiple cigars within the bed region we continue with our loop
			if s.Cigar[j].Op == cigar.Match {
				ans = append(ans, s.Seq[idxStart:idx]...)
			} else if s.Cigar[j].Op == cigar.Insertion {
				tmp = s.Seq[idxStart:idx]
				dna.AllToLower(tmp)
				ans = append(ans, tmp...)
			} else if cigar.ConsumesReference(s.Cigar[j].Op) {
				ans = addGaps(ans, idx-idxStart)
			}
		}
	}

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

func addGaps(ans []dna.Base, numGap int) []dna.Base {
	for i := 0; i < numGap; i++ {
		ans = append(ans, dna.Gap)
	}
	return ans
}
