package maf

import (
	"sort"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
)

// MafSLinesAreEqual returns true if, and only if, the values
// in two "S" lines are identical.
func MafSLinesAreEqual(a MafSLine, b MafSLine) bool {
	if a.Src != b.Src {
		return false
	}
	if a.Start != b.Start {
		return false
	}
	if a.Size != b.Size {
		return false
	}
	if a.Strand != b.Strand {
		return false
	}
	if a.SrcSize != b.SrcSize {
		return false
	}
	if dna.CompareSeqsCaseSensitive(a.Seq, b.Seq) != 0 {
		return false
	}
	return true
}

// MafILinesAreEqual returns true if, and only if, the values
// in two "I" lines are identical.
func MafILinesAreEqual(a MafILine, b MafILine) bool {
	if a.Src != b.Src {
		return false
	}
	if a.LeftStatus != b.LeftStatus {
		return false
	}
	if a.LeftCount != b.LeftCount {
		return false
	}
	if a.RightStatus != b.RightStatus {
		return false
	}
	if a.RightCount != b.RightCount {
		return false
	}
	return true
}

// MafELinesAreEqual returns true if, and only if, the values
// in two "E" lines are identical.
func MafELinesAreEqual(a MafELine, b MafELine) bool {
	if a.Src != b.Src {
		return false
	}
	if a.Start != b.Start {
		return false
	}
	if a.Size != b.Size {
		return false
	}
	if a.Strand != b.Strand {
		return false
	}
	if a.SrcSize != b.SrcSize {
		return false
	}
	if a.Status != b.Status {
		return false
	}
	return true
}

// MafSpeciesAreEqual returns true if, and only if, the
// info for two MafSpecies structs are identical.  This
// is the Src field, as well as any "S," "I," or "E" lines.
func MafSpeciesAreEqual(a MafSpecies, b MafSpecies) bool {
	if a.Src != b.Src {
		return false
	}
	if !MafSLinesAreEqual(*a.SLine, *b.SLine) { //a.SLine has type *MafSLine, so compare &a.SLine
		return false
	}
	if !MafILinesAreEqual(*a.ILine, *b.ILine) {
		return false
	}
	if !MafELinesAreEqual(*a.ELine, *b.ELine) {
		return false
	}
	return true
}

// MafsAreEqual returns true if, and only if, two maf
// blocks have identical values, including the species
// in each block and any "S," "I," or "E" lines associated
// with each species.
func MafsAreEqual(a Maf, b Maf) bool {
	if a.Score != b.Score {
		return false
	}
	if len(a.Species) != len(b.Species) {
		return false
	}
	for i := 0; i < len(a.Species); i++ {
		if !MafSpeciesAreEqual(*a.Species[i], *b.Species[i]) {
			return false
		}
	}
	return true
}

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

// SortByPos will sort a slice of maf blocks by their genomic position,
// which will first use the value of Src (i.e. assembly and/or chromosome)
// and then the numerical position.
func SortByPos(m []*Maf) {
	sort.Slice(m, func(i, j int) bool { return comparePos(m[i], m[j]) == -1 })
}

// SortByPos will sort a slice of maf blocks in reverse order (highest to lowest) by their genomic position,
// which will first use the value of Src (i.e. assembly and/or chromosome)
// and then the numerical position.
func SortByPosRev(m []*Maf) {
	sort.Slice(m, func(i, j int) bool { return comparePos(m[i], m[j]) == 1 })
}

// SortScore will sort a slice of maf blocks by their score.
func SortScore(m []*Maf) {
	sort.Slice(m, func(i, j int) bool { return compareScore(m[i], m[j]) == 1 })
}
