package psl

import (
	"github.com/vertgenlab/gonomics/numbers"
)

// Equal will check two psl structs and compare every single field of data to determine if they are the same.
func Equal(x, y Psl) bool {
	if x.Match != y.Match {
		return false
	}
	if x.MisMatch != y.MisMatch {
		return false
	}
	if x.RepeatMatch != y.RepeatMatch {
		return false
	}
	if x.Ns != y.Ns {
		return false
	}
	if x.QNumIns != y.QNumIns {
		return false
	}
	if x.QBaseIns != y.QBaseIns {
		return false
	}
	if x.TNumIns != y.TNumIns {
		return false
	}
	if x.TBaseIns != y.TBaseIns {
		return false
	}
	if x.Strand != y.Strand {
		return false
	}
	if x.QName != y.QName {
		return false
	}
	if x.QSize != y.QSize {
		return false
	}
	if x.QStart != y.QStart {
		return false
	}
	if x.QEnd != y.QEnd {
		return false
	}
	if x.TName != y.TName {
		return false
	}
	if x.TSize != y.TSize {
		return false
	}
	if x.TStart != y.TStart {
		return false
	}
	if x.BlockCount != y.BlockCount {
		return false
	}
	if !numbers.EqualSliceInt(x.BlockSize, y.BlockSize) {
		return false
	}
	if !numbers.EqualSliceInt(x.QList, y.QList) {
		return false
	}
	if !numbers.EqualSliceInt(x.TList, y.TList) {
		return false
	}
	return true
}
