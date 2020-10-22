package genePred

import (
	"log"
	"strings"
)

func AllAreEqual(a []*GenePred, b []*GenePred) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if !Equal(a[i], b[i]) {
			return false
		}
	}
	return true
}

func Equal(a *GenePred, b *GenePred) bool {
	if strings.Compare(a.Id, b.Id) != 0 {
		log.Print("1")
		return false
	}
	if strings.Compare(a.Chrom, b.Chrom) != 0 {
		log.Print("2")
		return false
	}
	if a.Strand != b.Strand {
		log.Print("3")
		return false
	}
	if a.TxStart != b.TxStart {
		log.Print("4")
		return false
	}
	if a.TxEnd != b.TxEnd {
		log.Print("5")
		return false
	}
	if a.CdsStart != b.CdsStart {
		log.Print("6")
		return false
	}
	if a.CdsEnd != b.CdsEnd {
		log.Print("7")
		return false
	}
	//exon ends must have the same number of values as exon starts
	for i := 0; i < len(a.ExonStarts); i++ {
		if a.ExonStarts[i] != b.ExonStarts[i] {
			log.Print("8")
			return false
		}
		if b.ExonEnds[i] != b.ExonEnds[i] {
			log.Print("9")
			return false
		}
	}
	return true
}
