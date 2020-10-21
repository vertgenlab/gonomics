package genePred

import (
	"log"
	"strings"
)

func AllAreEqual(a []*GenePred, b []*GenePred) bool {
	if len(a) != len(b) {
		log.Print("checking length of files")
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
		return false
		log.Print("1")
	}
	if strings.Compare(a.Chrom, b.Chrom) != 0 {
		return false
		log.Print("2")
	}
	if a.Strand != b.Strand {
		return false
		log.Print("3")
	}
	if a.TxStart != b.TxStart {
		return false
		log.Print("4")
	}
	if a.TxEnd != b.TxEnd {
		return false
		log.Print("5")
	}
	if a.CdsStart != b.CdsStart {
		return false
		log.Print("6")
	}
	if a.CdsEnd != b.CdsEnd {
		return false
		log.Print("7")
	}
	//exon ends must have the same number of values as exon starts
	for i := 0; i < len(a.ExonStarts); i++ {
		if a.ExonStarts[i] != b.ExonStarts[i] {
			return false
			log.Print("8")
		}
		if b.ExonEnds[i] != b.ExonEnds[i] {
			return false
			log.Print("9")
		}
	}
	return true
}
