package chain

import (
	"log"
)

// TPosToQPos converts a target position in a chain to the corresponding query position.
// The bool returns true if the returned TPos was found in a (Size) section of a Chain, and
// false otherwise, as in when the TPos fell in TBases.
func TPosToQPos(c Chain, TPos int) (int, bool) {
	var currT, currQ int

	//if target is negative strand, we convert TStart to the reverse complement position.
	if c.TStrand {
		currT = c.TStart
	} else {
		log.Fatalf("All target strands should be positive. Current chain:\n %s\n", ToString(c))
	}

	if c.QStrand {
		currQ = c.QStart
	} else {
		currQ = c.QEnd - 1
	}

	if TPos < c.TStart || TPos > c.TEnd {
		log.Fatalf("Error in TPosToQPos: TPos: %d, is not within the range of the chain. TStart: %d. TEnd: %d.", TPos, c.TStart, c.TEnd)
	}
	for i := 0; i < len(c.Alignment); i++ {
		if c.QStrand {
			if currT+c.Alignment[i].Size > TPos {
				return currQ + (TPos - currT), true
			}
			currT += c.Alignment[i].Size
			currQ += c.Alignment[i].Size
			if currT+c.Alignment[i].TBases > TPos {
				//in this case, the TPos is in a place with no corresponding query block.
				return currQ, false
			}
			currT += c.Alignment[i].TBases
			currQ += c.Alignment[i].QBases
		} else {
			if currT + c.Alignment[i].Size > TPos {
				return currQ - (TPos - currT), true
			}
			currT += c.Alignment[i].Size
			currQ -= c.Alignment[i].Size
			if currT+c.Alignment[i].Size > TPos {
				return currQ, false
			}
		}
	}

	log.Fatalf("Unable to locate the TPos within chain.")
	return -1, false
}
