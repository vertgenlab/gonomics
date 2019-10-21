package simpleGraph

import (
//	"fmt"
//	"github.com/vertgenlab/gonomics/common"
//	"github.com/vertgenlab/gonomics/dna"
//	"github.com/vertgenlab/gonomics/fileio"
//	"io"
//	"strings"
)

type Seed struct {
	Id  int64
	Start int64
	End int64
	Score int64
}
/*
func overlapSeed(a *Seed, newId, newStart, newEnd) bool {
	if a.Id == newId && common.MaxInt64(a.Start, newStart) < common.MinInt64(a.End, newEnd) {
		return true
	} else {
		return false
	}
}

func addSeed(existing []*Seed, newId, newStart, newEnd) []*Seed {
	for i:=0; i < len(existing); i++ {
		if overlapSeed(existing[i], newId, newStart, newEnd) {
			existing[i].Start = common.MinInt64(existing[i].Start, newStart)
			existing[i].End = common.MaxInt64(existing[i].End, newEnd)
			existing[i].Score++
		}
}
*/
