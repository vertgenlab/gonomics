package axt

import (
	"sort"
	"strings"
)

func compareRCoord(alpha Axt, beta Axt) int {
	if alpha.RStart < beta.RStart {
		return -1
	}
	if alpha.RStart > beta.RStart {
		return 1
	}
	if alpha.REnd < beta.REnd {
		return -1
	}
	if alpha.REnd > beta.REnd {
		return 1
	}
	return 0
}

func compareRName(alpha Axt, beta Axt) int {
	return strings.Compare(alpha.RName, beta.RName)
}

func compareQName(alpha Axt, beta Axt) int {
	return strings.Compare(alpha.QName, beta.QName)
}

func compareScore(alpha Axt, beta Axt) int {
	if alpha.Score < beta.Score {
		return 1
	}
	if alpha.Score > beta.Score {
		return -1
	}
	return 0
}

func compareRNameCoord(alpha Axt, beta Axt) int {
	compareStorage := compareRName(alpha, beta)
	if compareStorage != 0 {
		return compareStorage
	} else {
		return compareRCoord(alpha, beta)
	}
}

// SortByRNameCoord sorts a slice of Axts by the reference chrom name and then genomic position
func SortByRNameCoord(axts []Axt) {
	sort.Slice(axts, func(i, j int) bool { return compareRNameCoord(axts[i], axts[j]) == -1 })
}

// SortByScore sorts a slice of Axts by their score.
func SortByScore(axts []Axt) {
	sort.Slice(axts, func(i, j int) bool { return compareScore(axts[i], axts[j]) == -1 })
}

func isEqual(alpha Axt, beta Axt) bool {
	if strings.Compare(alpha.RName, beta.RName) != 0 {
		return false
	}
	if alpha.RStart != beta.RStart {
		return false
	}
	if alpha.REnd != beta.REnd {
		return false
	}
	if strings.Compare(alpha.QName, beta.QName) != 0 {
		return false
	}
	if alpha.QStart != beta.QStart {
		return false
	}
	if alpha.QEnd != beta.QEnd {
		return false
	}
	if alpha.QStrandPos != beta.QStrandPos {
		return false
	}
	if alpha.Score != beta.Score {
		return false
	}
	if len(alpha.RSeq) != len(beta.RSeq) {
		return false
	}
	for i := 0; i < len(alpha.RSeq); i++ {
		if alpha.RSeq[i] != beta.RSeq[i] {
			return false
		}
	}
	if len(alpha.QSeq) != len(beta.QSeq) {
		return false
	}
	for j := 0; j < len(alpha.QSeq); j++ {
		if alpha.QSeq[j] != beta.QSeq[j] {
			return false
		}
	}
	return true
}

func allEqual(alpha []Axt, beta []Axt) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := 0; i < len(alpha); i++ {
		if !isEqual(alpha[i], beta[i]) {
			return false
		}
	}
	return true
}
