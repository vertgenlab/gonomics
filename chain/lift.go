package chain

import (
	"log"
)

//TPosToQPos converts a target position in a chain to the corresponding query position.
func TPosToQPos(c Chain, TPos int) QPos {
	if TPos < c.TStart || TPos > c.TEnd {
		log.Fatalf("Error in TPosToQPos: TPos: %d, is not within the range of the chain. TStart: %d. TEnd: %d.", TPos, c.TStart, c.TEnd)
	}
	var currT int = c.TStart
	var currQ int = c.QStart

	for i := 0; i < len(c.Alignments); i++ {
		if currT + c.Alignments[i].Size > TPos {
			return currQ + TPos - currT
		}
		currT += c.Alignments[i].Size
		currQ += c.Alignments[i].Size
		if currT + c.Alignments[i].TBases > TPos {
			//in this case, the TPos is in a place with no corresponding query block.
			//Should we just return currQ?
			return currQ
		}
		currT += c.Alignments[i].TBases
		currQ += c.Alignments[i].QBases
	}
	log.Fatalf("Unable to locate the TPos within chain.")
	return -1
}

