package maf

import (
	"sort"
	"strings"
)

func comparePos(a *Maf, b *Maf) int {
	var aSLine, bSLine *MafSLine
	aSLine = a.Species[0].SLine
	bSLine = b.Species[0].SLine
        chromComp := strings.Compare(aSLine.Src, bSLine.Src)
        if chromComp != 0 {
                return chromComp
        } else if aSLine.Start < bSLine.Start {
                return -1
        } else if aSLine.Start > bSLine.Start {
                return 1
        } else if aSLine.Size < bSLine.Size {
                return -1
        } else if aSLine.Size > bSLine.Size {
                return 1
        } else {
		return 0
	}
}

func compareScore(a *Maf, b *Maf) int {
        if a.Score < b.Score {
                return -1
        } else if a.Score > b.Score {
                return 1
        } else {
		return 0
	}
}

func SortByPos(m []*Maf) {
        sort.Slice(m, func(i, j int) bool { return comparePos(m[i], m[j]) == -1 })
}

func SortByPosRev(m []*Maf) {
	sort.Slice(m, func(i, j int) bool { return comparePos(m[i], m[j]) == 1 })
}

func SortScore(m []*Maf) {
        sort.Slice(m, func(i, j int) bool { return compareScore(m[i], m[j]) == 1 })
}

