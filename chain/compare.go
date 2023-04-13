package chain

import (
	"log"
	"sort"
	"strings"
)

// Uses bool to compare target or query coordinates as one function
// true is for target, false is for query
func compareStartCoord(a Chain, b Chain, checkTarget bool) int {
	if checkTarget {
		return compareTargetCoord(a, b)
	} else {
		return compareQueryCoord(a, b)
	}
}

func compareTargetCoord(a Chain, b Chain) int {
	sameChr := strings.Compare(a.TName, b.TName)
	if sameChr != 0 {
		return sameChr
	}
	var alphaStart, betaStart int
	//check for negative strand and swap if you have to
	if !a.TStrand {
		alphaStart = getSwapTCoord(a, true, false)
	} else {
		alphaStart = a.TStart
	}
	//check for negative again
	if !b.TStrand {
		betaStart = getSwapTCoord(b, true, false)
	} else {
		betaStart = b.TStart
	}
	//finally, we do the comparison
	if alphaStart < betaStart {
		return -1
	}
	if alphaStart > betaStart {
		return 1
	}
	return 0
}

func compareQueryCoord(a Chain, b Chain) int {
	sameChr := strings.Compare(a.QName, b.QName)
	if sameChr != 0 {
		return sameChr
	}
	var alphaStart, betaStart int
	//check for negative strand and swap if you have to
	if !a.QStrand {
		alphaStart = getSwapQCoord(a, true, false)
	} else {
		alphaStart = a.QStart
	}
	//check for negative again
	if !b.QStrand {
		betaStart = getSwapQCoord(b, true, false)
	} else {
		betaStart = b.QStart
	}

	if alphaStart < betaStart {
		return -1
	}
	if alphaStart > betaStart {
		return 1
	}

	return 0
}

// true/false bool to either sort by target or query
// true=target, false=query
func SortByCoordinates(align []Chain, whichGenome bool) {
	sort.Slice(align, func(i, j int) bool { return compareStartCoord(align[i], align[j], whichGenome) == -1 })
}

func compareScores(a Chain, b Chain) int {
	if a.Score < b.Score {
		return -1
	}
	if a.Score > b.Score {
		return 1
	}
	return 0
}

func Equal(a []Chain, b []Chain) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !isEqual(a[i], b[i]) {
			return false
		}
	}
	return true
}

func isEqual(a Chain, b Chain) bool {
	if strings.Compare(a.TName, b.TName) != 0 {
		return false
	}
	if strings.Compare(a.QName, b.QName) != 0 {
		return false
	}
	if a.TSize != b.TSize {
		return false
	}
	if a.TStrand != b.TStrand {
		return false
	}
	if a.TStart != b.TStart {
		return false
	}
	if a.TEnd != b.TEnd {
		return false
	}
	if a.QSize != b.QSize {
		return false
	}
	if a.QStrand != b.QStrand {
		return false
	}
	if a.QStart != b.QStart {
		return false
	}
	if a.QEnd != b.QEnd {
		return false
	}
	if a.Id != b.Id {
		return false
	}
	return true
}

// set start true if you want to adjust the start or end to be true to adjust for the negative stand on the end
func getSwapTCoord(ch Chain, start bool, end bool) int {
	if !ch.TStrand {
		if start {
			return ch.TSize - ch.TEnd
		}
		if end {
			return ch.TSize - ch.TStart
		}
		if start && end {
			log.Fatalf("Error: Must select either start or end, and not both...\n")
		}
	} else {
		log.Fatalf("Error: Positive strand chain record detected, This function is primarily used only to swap coordinates on the negative strand ...\n")
	}
	return -1
}

// Could have one function that performs on both query and target
func getSwapQCoord(ch Chain, start bool, end bool) int {
	if !ch.QStrand {
		if start {
			return ch.QSize - ch.QEnd
		}
		if end {
			return ch.QSize - ch.QStart
		}
		if start && end {
			log.Fatalf("Error: Must select either start or end, and not both...\n")
		}
	} else {
		log.Fatalf("Error: Positive strand chain record detected, This function is primarily used only to swap coordinates on the negative strand ...\n")

	}
	return -1
}
