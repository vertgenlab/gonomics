package chain

import (
	"log"
)

//TPosToQPos converts a target position in a chain to the corresponding query position.
func TPosToQPos(c *Chain, TPos int) int {
	if TPos < c.TStart || TPos > c.TEnd {
		log.Fatalf("Error in TPosToQPos: TPos: %d, is not within the range of the chain. TStart: %d. TEnd: %d.", TPos, c.TStart, c.TEnd)
	}
	var currT int = c.TStart
	var currQ int = c.QStart

	for i := 0; i < len(c.Alignment); i++ {
		if c.TStrand {
			if currT + c.Alignment[i].Size > TPos {
				if c.QStrand {
					return currQ + TPos - currT
				} else {
					return currQ - (TPos - currT)
				}
			}
			currT += c.Alignment[i].Size
			if c.QStrand {
				currQ += c.Alignment[i].Size
			} else {
				currQ -= c.Alignment[i].Size
			}
		} else {

		}




		currT += c.Alignment[i].Size
		currQ += c.Alignment[i].Size
		if currT + c.Alignment[i].TBases > TPos {
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


/* Old
func TPosToQPos(c Chain, TPos int) int {
	if TPos < c.TStart || TPos > c.TEnd {
		log.Fatalf("Error in TPosToQPos: TPos: %d, is not within the range of the chain. TStart: %d. TEnd: %d.", TPos, c.TStart, c.TEnd)
	}
	var currT int = c.TStart
	var currQ int = c.QStart

	for i := 0; i < len(c.Alignment); i++ {
		if currT + c.Alignment[i].Size > TPos {
			return currQ + TPos - currT
		}
		currT += c.Alignment[i].Size
		currQ += c.Alignment[i].Size
		if currT + c.Alignment[i].TBases > TPos {
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
*/
