package chain

import (
	"log"
)

//TPosToQPos converts a target position in a chain to the corresponding query position.
func TPosToQPos(c *Chain, TPos int) int {
	var currT, currQ int

	//if target is negative strand, we convert TStart to the reverse complement position.
	if c.TStrand {
		currT = c.TStart
	} else {
		currT = c.TSize - 1 - c.TStart
	}

	if c.QStrand {
		currQ = c.QStart
	} else {
		currQ = c.QSize - 1 - c.QStart
	}

	if TPos < c.TStart || TPos > c.TEnd {
		log.Fatalf("Error in TPosToQPos: TPos: %d, is not within the range of the chain. TStart: %d. TEnd: %d.", TPos, c.TStart, c.TEnd)
	}

	for i := 0; i < len(c.Alignment); i++ {
		if currT+c.Alignment[i].Size > TPos {
			return currQ + TPos - currT
		}
		currT += c.Alignment[i].Size
		currQ += c.Alignment[i].Size
		if currT+c.Alignment[i].TBases > TPos {
			//in this case, the TPos is in a place with no corresponding query block.
			//Should we just return currQ?
			return currQ
		}
		currT += c.Alignment[i].TBases
		currQ += c.Alignment[i].QBases
	}

	log.Fatalf("Unable to locate the TPos within chain.")
	return -1
}
