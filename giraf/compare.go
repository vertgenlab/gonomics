package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"strings"
)

func isEqual(a *Giraf, b *Giraf) bool {
	if strings.Compare(a.QName, b.QName) != 0 {
		return false
	}
	if a.QStart != b.QStart {
		return false
	}
	if a.QEnd != b.QEnd {
		return false
	}
	if a.PosStrand != b.PosStrand {
		return false
	}
	if strings.Compare(PathToString(a.Path), PathToString(b.Path)) != 0 {
		return false
	}
	if strings.Compare(cigar.ToString(a.Aln), cigar.ToString(b.Aln)) != 0 {
		return false
	}
	if a.AlnScore != b.AlnScore {
		return false
	}
	if a.MapQ != b.MapQ {
		return false
	}
	if dna.CompareSeqsIgnoreCase(a.Seq, b.Seq) != 0 {
		return false
	}
	//TODO: write a check for qual
	//if a.Qual != b.Qual
	if strings.Compare(NotesToString(a.Notes), NotesToString(b.Notes)) != 0 {
		return false
	}
	return true
}

func AllEqual(a []*Giraf, b []*Giraf) bool {
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
