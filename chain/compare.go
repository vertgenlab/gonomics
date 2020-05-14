package chain

import (
	"strings"
)

func isEqual(a *Chain, b *Chain) bool {
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
