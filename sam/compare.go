package sam

import (
	"sort"
	"strings"
)

func compareName(alpha *SamAln, beta *SamAln) int {
	return strings.Compare(alpha.RName, beta.RName)
}

func SortByName(alignments []*SamAln) {
	sort.Slice(alignments, func(i, j int) bool { return compareName(alignments[i], alignments[j]) == -1 })
}

func SortByCoord(samRecords []*SamAln) {
	sort.Slice(samRecords, func(i, j int) bool { return Compare(samRecords[i], samRecords[j]) == -1 })
}

func Compare(a *SamAln, b *SamAln) int {
	chrName := strings.Compare(a.RName, b.RName)
	if chrName != 0 {
		return chrName
	}
	if a.Pos < b.Pos {
		return -1
	}
	if a.Pos > b.Pos {
		return 1
	}
	return 0
}

func SplitSamByChr(samRecords []*SamAln) map[string][]*SamAln {
	SortByCoord(samRecords)
	genome := make(map[string][]*SamAln)
	for _, alignedRead := range samRecords {
		genome[alignedRead.RName] = append(genome[alignedRead.RName], alignedRead)
	}
	return genome
}
