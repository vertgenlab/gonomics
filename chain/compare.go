package chain

import (
	"sort"
	"strings"
)

//Uses bool to compare target or query coordinates as one function
//true is for target, false is for query
func compareStartCoord(a *Chain, b *Chain, whichGenome bool) int {
	if whichGenome {
		sameChr := strings.Compare(a.TName, b.TName)
		if sameChr != 0 {
			return sameChr
		}
		if a.TStart < b.TStart {
			return -1
		}
		if a.TStart > b.TStart {
			return 1
		}
	} else {
		sameChr := strings.Compare(a.QName, b.QName)
		if sameChr != 0 {
			return sameChr
		}
		if a.QStart < b.QStart {
			return -1
		}
		if a.QStart > b.QStart {
			return 1
		}
	}
	return 0
}

//true/false bool to either sort by target or query
//true=target, false=query
func SortByCoordinates(align []*Chain, whichGenome bool) {
	sort.Slice(align, func(i, j int) bool { return compareStartCoord(align[i], align[j], whichGenome) == -1 })
}

func compareScores(a *Chain, b *Chain) int {
	if a.Score < b.Score {
		return -1
	}
	if a.Score > b.Score {
		return 1
	}
	return 0
}

func CompareName(alpha string, beta string) int {
	return strings.Compare(alpha, beta)
}

func Equal(a []*Chain, b []*Chain) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !sameRecord(a[i], b[i]) {
			return false
		}
	}
	return true
}

func sameRecord(a *Chain, b *Chain) bool {
	if strings.Compare(a.TName, b.TName) != 0 {
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
	if strings.Compare(a.QName, b.QName) != 0 {
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
